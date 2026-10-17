<#
.SYNOPSIS
Times bench_odr built with reference BLAS and with OpenBLAS.

.DESCRIPTION
Requires builddir-refblas and builddir-openblas configured with -Dbuild_benchmark=true.
Writes raw timings to benchmark\results\raw.csv and per-case medians to benchmark\results\summary.csv.
#>
param(
    [int]$Reps = 10,
    [string]$OpenblasBin = 'C:\Users\hugom\micromamba\envs\fortran\Library\bin',
    [string]$OutDir = 'benchmark\results'
)

$ErrorActionPreference = 'Stop'
Set-Location (Join-Path $PSScriptRoot '..')
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
$env:Path = "$OpenblasBin;$env:Path"

$configs = @(
    [pscustomobject]@{ Label = 'refblas';     Exe = 'builddir-refblas\benchmark\bench_odr.exe';  Threads = $null }
    [pscustomobject]@{ Label = 'openblas_mt'; Exe = 'builddir-openblas\benchmark\bench_odr.exe'; Threads = $null }
    [pscustomobject]@{ Label = 'openblas_1t'; Exe = 'builddir-openblas\benchmark\bench_odr.exe'; Threads = '1' }
)

function Invoke-Bench {
    param($Config, [int]$Rep)
    if ($Config.Threads) { $env:OPENBLAS_NUM_THREADS = $Config.Threads }
    else { Remove-Item Env:OPENBLAS_NUM_THREADS -ErrorAction SilentlyContinue }
    $lines = & $Config.Exe $Rep
    if ($LASTEXITCODE -ne 0) { throw "bench_odr failed for $($Config.Label) (exit $LASTEXITCODE)" }
    $lines | Where-Object { $_ -match '^(linear|expo),' }
}

function Get-Median {
    param([double[]]$Values)
    $s = @($Values | Sort-Object)
    $k = $s.Count
    if ($k % 2) { return $s[[int](($k - 1) / 2)] }
    return ($s[[int]($k / 2) - 1] + $s[[int]($k / 2)]) / 2
}

# Warm-up run (not recorded)
foreach ($c in $configs) { Invoke-Bench $c 0 | Out-Null }

$rows = [System.Collections.Generic.List[object]]::new()
for ($rep = 1; $rep -le $Reps; $rep++) {
    # Alternate run order so that drift affects each build equally
    $order = if ($rep % 2 -eq 1) { $configs } else { $configs[($configs.Count - 1)..0] }
    foreach ($c in $order) {
        foreach ($line in (Invoke-Bench $c $rep)) {
            $f = $line.Split(',')
            $rows.Add([pscustomobject]@{
                Config = $c.Label
                Model  = $f[0]
                N      = [int]$f[1]
                Bound  = $f[2].Trim()
                Rep    = [int]$f[3]
                Time_s = [double]($f[4].Trim())
                Info   = [int]$f[5]
                Beta   = $f[6].Trim()
            })
        }
    }
    Write-Host "Repetition $rep/$Reps done"
}

$rows | Export-Csv -NoTypeInformation -Path (Join-Path $OutDir 'raw.csv')

$cases = $rows | Where-Object { $_.Config -eq 'refblas' } | Select-Object Model, N, Bound -Unique
$summary = foreach ($cs in $cases) {
    $key = { param($cfg) $rows | Where-Object { $_.Config -eq $cfg -and $_.Model -eq $cs.Model -and $_.N -eq $cs.N -and $_.Bound -eq $cs.Bound } }
    $med = @{}
    foreach ($c in $configs) {
        $med[$c.Label] = Get-Median ([double[]]@((& $key $c.Label) | ForEach-Object { $_.Time_s }))
    }
    $ref1 = (& $key 'refblas') | Where-Object { $_.Rep -eq 1 } | Select-Object -First 1
    $relMT = 0.0
    $rel1T = 0.0
    $b0 = [double[]]($ref1.Beta -split ';')
    foreach ($label in 'openblas_mt', 'openblas_1t') {
        $b1 = [double[]](((& $key $label) | Where-Object { $_.Rep -eq 1 } | Select-Object -First 1).Beta -split ';')
        $maxRel = 0.0
        for ($i = 0; $i -lt $b0.Count; $i++) {
            $d = [Math]::Abs($b0[$i] - $b1[$i]) / [Math]::Max([Math]::Abs($b0[$i]), 1e-300)
            if ($d -gt $maxRel) { $maxRel = $d }
        }
        if ($label -eq 'openblas_mt') { $relMT = $maxRel } else { $rel1T = $maxRel }
    }
    [pscustomobject]@{
        Model          = $cs.Model
        N              = $cs.N
        Bound          = $cs.Bound
        Info           = $ref1.Info
        RefBLAS_s      = [Math]::Round($med['refblas'], 6)
        OpenBLAS_MT_s  = [Math]::Round($med['openblas_mt'], 6)
        OpenBLAS_1T_s  = [Math]::Round($med['openblas_1t'], 6)
        Saving_MT_pct  = [Math]::Round(100 * (1 - $med['openblas_mt'] / $med['refblas']), 1)
        Saving_1T_pct  = [Math]::Round(100 * (1 - $med['openblas_1t'] / $med['refblas']), 1)
        MaxRelBetaDiff_MT = '{0:E2}' -f $relMT
        MaxRelBetaDiff_1T = '{0:E2}' -f $rel1T
    }
}

$totals = [pscustomobject]@{
    Model          = 'TOTAL'
    N              = ''
    Bound          = ''
    Info           = ''
    RefBLAS_s      = [Math]::Round(($summary | Measure-Object RefBLAS_s -Sum).Sum, 6)
    OpenBLAS_MT_s  = [Math]::Round(($summary | Measure-Object OpenBLAS_MT_s -Sum).Sum, 6)
    OpenBLAS_1T_s  = [Math]::Round(($summary | Measure-Object OpenBLAS_1T_s -Sum).Sum, 6)
    Saving_MT_pct  = [Math]::Round(100 * (1 - ($summary | Measure-Object OpenBLAS_MT_s -Sum).Sum / ($summary | Measure-Object RefBLAS_s -Sum).Sum), 1)
    Saving_1T_pct  = [Math]::Round(100 * (1 - ($summary | Measure-Object OpenBLAS_1T_s -Sum).Sum / ($summary | Measure-Object RefBLAS_s -Sum).Sum), 1)
    MaxRelBetaDiff_MT = ''
    MaxRelBetaDiff_1T = ''
}
$summary = @($summary) + $totals

$summary | Export-Csv -NoTypeInformation -Path (Join-Path $OutDir 'summary.csv')
$summary | Format-Table -AutoSize
