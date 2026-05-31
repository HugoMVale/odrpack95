#include <stdio.h>

#include "../include/odrpack/odrpack.h"

int main() {
    int info_values[] = {1, 2, 3, 4, 20, 21, 300, 4000, 40000, 50000, 60002, 70000, 80100,
                         90100, 32132132};
    int nvalues = (int)(sizeof(info_values) / sizeof(info_values[0]));

    for (int i = 0; i < nvalues; i++) {
        char message[256] = {0};
        int info = info_values[i];
        stop_message_c(info, (int)sizeof(message), message);
        printf("Stop reason (info = %d): %s\n", info, message);
    }

    return 0;
}
