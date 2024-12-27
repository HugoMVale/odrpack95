#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "odrpack.h"

namespace py = pybind11;
using namespace pybind11::literals;

// // Wrapper for the odr_short_c function
// void odr_short_c_wrapper(
//     py::object fcn,
//     int n,
//     int m,
//     int np,
//     int nq,
//     py::array_t<double> beta,
//     py::array_t<double> y,
//     py::array_t<double> x,
//     py::array_t<double> delta,
//     py::array_t<double> lower,
//     py::array_t<double> upper,
//     int job)
// {
//     auto beta_ptr = beta.mutable_data();
//     auto y_ptr = y.data();
//     auto x_ptr = x.data();
//     auto delta_ptr = delta.mutable_data();
//     auto lower_ptr = lower.data();
//     auto upper_ptr = upper.data();

//     // Define a C++ lambda as the callback
//     auto c_fcn = [fcn](int *n, int *m, int *np, int *nq, int *ldn, int *ldm, int *ldnp,
//                        double *beta, double *xplusd, int *ifixb, int *ifixx, int *ldifx, int *ideval,
//                        double *f, double *fjacb, double *fjacd, int *istop)
//     {
//         py::gil_scoped_acquire gil; // Acquire the GIL since we are calling Python
//         py::array_t<double> beta_py({*np}, beta);
//         py::array_t<double> xplusd_py({*m, *ldn}, xplusd);
//         py::array_t<int> ifixb_py({*np}, ifixb);
//         py::array_t<int> ifixx_py({*m, *ldifx}, ifixx);
//         py::array_t<double> f_py({*nq, *ldn}, f);
//         py::array_t<double> fjacb_py({*nq, *ldnp, *ldn}, fjacb);
//         py::array_t<double> fjacd_py({*nq, *ldm, *ldn}, fjacd);

//         py::object ret = fcn(*n, *m, *np, *nq, beta_py, xplusd_py, ifixb_py, ifixx_py, *ideval, f_py, fjacb_py, fjacd_py);
//         *istop = ret.cast<int>();
//     };

//     // Call the C function
//     odr_short_c(c_fcn, &n, &m, &np, &nq, beta_ptr, y_ptr, x_ptr, delta_ptr, lower_ptr, upper_ptr, &job);
// }

// PYBIND11_MODULE(_odrpack, m)
// {
//     m.def("odr_short_c", &odr_short_c_wrapper,
//           "Wrapper for the odr_short_c function.",
//           py::arg("fcn"),
//           py::arg("n"),
//           py::arg("m"),
//           py::arg("np"),
//           py::arg("nq"),
//           py::arg("beta"),
//           py::arg("y"),
//           py::arg("x"),
//           py::arg("delta"),
//           py::arg("lower"),
//           py::arg("upper"),
//           py::arg("job"), );
// }

PYBIND11_MODULE(_odrpack, m)
{
  m.def("open_file", [](int lun, const char *filename)
        {
        int ierr = 0;
        open_file(&lun, filename, &ierr);
        return std::make_tuple(lun, ierr); },

        "Open a new file associated with a specified logical unit number.\n"
        "\n"
        "Arguments:\n"
        "    lun      (int): Logical unit number.\n"
        "    filename (str): String containing the file name.\n"
        "Returns:\n"
        "    tuple: A tuple containing the updated lun (int) and ierr (int).",

        py::arg("lun"), py::arg("filename"));

  m.def("close_file", [](int lun)
        {
        int ierr = 0;
        close_file(&lun, &ierr);
        return ierr; },

        "Close a file associated with a specified logical unit number.\n"
        "\n"
        "Arguments:\n"
        "    lun    (int): Logical unit number.\n"
        "Returns:\n"
        "    ierr   (int): Error code (compiler dependent).",

        py::arg("lun"));

  m.def("workspace_dimensions", [](int n, int m, int np, int nq, bool isodr)
        {
        int lwork = 0;
        int liwork = 0;
        workspace_dimensions_c(&n, &m, &np, &nq, &isodr, &lwork, &liwork);
        return std::make_tuple(lwork, liwork); },

        "Calculate the dimensions of the workspace arrays.\n"
        "\n"
        "Arguments:\n"
        "    n      (int): Number of observations.\n"
        "    m      (int): Number of columns of data in the explanatory variable.\n"
        "    np     (int): Number of function parameters.\n"
        "    nq     (int): Number of responses per observation.\n"
        "    isodr  (bool): Variable designating whether the solution is by ODR (True) or by OLS (False).\n"
        "Returns:\n"
        "    tuple: A tuple containing the lengths of the work arrays (lwork, liwork).",
        py::arg("n"), py::arg("m"), py::arg("np"), py::arg("nq"), py::arg("isodr"));

  m.def("diwinf", [](int n, int np, int nq)
        {
        iworkidx_t iworkidx = {};
        diwinf_c(&n, &np, &nq, &iworkidx);
        return py::dict(
            "msgb"_a = iworkidx.msgb,
            "msgd"_a = iworkidx.msgd,
            "ifix2"_a = iworkidx.ifix2,
            "istop"_a = iworkidx.istop,
            "nnzw"_a = iworkidx.nnzw,
            "npp"_a = iworkidx.npp,
            "idf"_a = iworkidx.idf,
            "job"_a = iworkidx.job,
            "iprin"_a = iworkidx.iprin,
            "luner"_a = iworkidx.luner,
            "lunrp"_a = iworkidx.lunrp,
            "nrow"_a = iworkidx.nrow,
            "ntol"_a = iworkidx.ntol,
            "neta"_a = iworkidx.neta,
            "maxit"_a = iworkidx.maxit,
            "niter"_a = iworkidx.niter,
            "nfev"_a = iworkidx.nfev,
            "njev"_a = iworkidx.njev,
            "int2"_a = iworkidx.int2,
            "irank"_a = iworkidx.irank,
            "ldtt"_a = iworkidx.ldtt,
            "bound"_a = iworkidx.bound,
            "liwkmn"_a = iworkidx.liwkmn
        ); },

        "Get storage locations within integer work space.\n"
        "\n"
        "Arguments:\n"
        "    n  (int): Number of observations.\n"
        "    np (int): Number of function parameters.\n"
        "    nq (int): Number of responses per observation.\n"
        "Returns:\n"
        "    dict: A dictionary containing the 0-based indexes of the integer work array.",

        py::arg("n"), py::arg("np"), py::arg("nq"));

  m.def("dwinf", [](int n, int m, int np, int nq, int ldwe, int ld2we, bool isodr)
        {
        workidx_t workidx = {};
        dwinf_c(&n, &m, &np, &nq, &ldwe, &ld2we, &isodr, &workidx);
        return py::dict(
            "delta"_a = workidx.delta,
            "eps"_a = workidx.eps,
            "xplus"_a = workidx.xplus,
            "fn"_a = workidx.fn,
            "sd"_a = workidx.sd,
            "vcv"_a = workidx.vcv,
            "rvar"_a = workidx.rvar,
            "wss"_a = workidx.wss,
            "wssde"_a = workidx.wssde,
            "wssep"_a = workidx.wssep,
            "rcond"_a = workidx.rcond,
            "eta"_a = workidx.eta,
            "olmav"_a = workidx.olmav,
            "tau"_a = workidx.tau,
            "alpha"_a = workidx.alpha,
            "actrs"_a = workidx.actrs,
            "pnorm"_a = workidx.pnorm,
            "rnors"_a = workidx.rnors,
            "prers"_a = workidx.prers,
            "partl"_a = workidx.partl,
            "sstol"_a = workidx.sstol,
            "taufc"_a = workidx.taufc,
            "epsma"_a = workidx.epsma,
            "beta0"_a = workidx.beta0,
            "betac"_a = workidx.betac,
            "betas"_a = workidx.betas,
            "betan"_a = workidx.betan,
            "s"_a = workidx.s,
            "ss"_a = workidx.ss,
            "ssf"_a = workidx.ssf,
            "qraux"_a = workidx.qraux,
            "u"_a = workidx.u,
            "fs"_a = workidx.fs,
            "fjacb"_a = workidx.fjacb,
            "we1"_a = workidx.we1,
            "diff"_a = workidx.diff,
            "delts"_a = workidx.delts,
            "deltn"_a = workidx.deltn,
            "t"_a = workidx.t,
            "tt"_a = workidx.tt,
            "omega"_a = workidx.omega,
            "fjacd"_a = workidx.fjacd,
            "wrk1"_a = workidx.wrk1,
            "wrk2"_a = workidx.wrk2,
            "wrk3"_a = workidx.wrk3,
            "wrk4"_a = workidx.wrk4,
            "wrk5"_a = workidx.wrk5,
            "wrk6"_a = workidx.wrk6,
            "wrk7"_a = workidx.wrk7,
            "lower"_a = workidx.lower,
            "upper"_a = workidx.upper,
            "lwkmn"_a = workidx.lwkmn
        ); },

        "Get storage locations within real work space.\n"
        "\n"
        "Arguments:\n"
        "    n (int): Number of observations.\n"
        "    m (int): Number of columns of data in the explanatory variable.\n"
        "    np (int): Number of function parameters.\n"
        "    nq (int): Number of responses per observation.\n"
        "    ldwe (int): Leading dimension of array `we`.\n"
        "    ld2we (int): Second dimension of array `we`.\n"
        "    isodr (bool): Solution by ODR (`True`) or OLS (`False`).\n"
        "Returns:\n"
        "    dict: A dictionary containing the 0-based indexes of the real work array.",

        py::arg("n"), py::arg("m"), py::arg("np"), py::arg("nq"), py::arg("ldwe"), py::arg("ld2we"), py::arg("isodr"));
}