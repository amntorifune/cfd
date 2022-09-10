#include <mpi.h>
#include <openacc.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

const int    N  = 1000;
const double h  = 1.0 / N;
const double k  = 1e3;
const double t0 = 0;
const double t1 = 100;

double *t    = nullptr;
double *tnew = nullptr;
double *t_global;

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int size, rank, ngpu;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ngpu = acc_get_num_devices(acc_device_nvidia);

    if (N % size != 0) {
        if (rank == 0) {
            printf("Err : N=%d is not divisible by -np=%d\n", N, size);
        }
        MPI_Finalize();
        return 0;
    }

    if (ngpu < size) {
        if (rank == 0) {
            printf("Err: ngpu=%d < size=%d\n", ngpu, size);
        }
        MPI_Finalize();
        return 0;
    }
    acc_set_device_num(rank, acc_device_nvidia);

    int len = N / size + 2;
    t    = new double[len];
    tnew = new double[len];
    memset(t   , 0, sizeof(double) * len);
    memset(tnew, 0, sizeof(double) * len);
    if (rank == 0) {
        t_global = new double[N];
        printf("running with %d procs\n", size);
    }

    double e;
    int    it = 0;
    MPI_Request req[4];

    clock_t start, end;
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        start = clock();
    }
    #pragma acc enter data copyin(t[0:len], tnew[0:len])
    do {
        double ee_sum_local  = 0;
        double ee_sum_global = 0;

        // exchange elements between subdomains
        if (size > 1) {
            if (rank == 0) {
                // only right end
                #pragma acc host_data use_device(t)
                MPI_Isend(&t[len - 2], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &req[0]);
                #pragma acc host_data use_device(t)
                MPI_Irecv(&t[len - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &req[1]);
            } else if (rank == size - 1) {
                // only left end
                #pragma acc host_data use_device(t)
                MPI_Isend(&t[1      ], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &req[2]);
                #pragma acc host_data use_device(t)
                MPI_Irecv(&t[0      ], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &req[3]);
            } else {
                // both right and left ends
                #pragma acc host_data use_device(t)
                MPI_Isend(&t[len - 2], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &req[0]);
                #pragma acc host_data use_device(t)
                MPI_Irecv(&t[len - 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &req[1]);
                #pragma acc host_data use_device(t)
                MPI_Isend(&t[1      ], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &req[2]);
                #pragma acc host_data use_device(t)
                MPI_Irecv(&t[0      ], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &req[3]);
            }
        }

        // inner elements, not envolving subdomains or boundaries
        #pragma acc kernels loop independent present(t[0:len], tnew[0:len]) copyin(len, k, h) reduction(+:ee_sum_local) copy(ee_sum_local)
        for (int i = 2; i <= len - 3; i ++) {
            double a_e, a_c, a_w, t_e, t_c, t_w, res, rhs;
            a_e = k / (h * h);
            a_c = - 2.0 * k / (h * h);
            a_w = k / (h * h);
            t_e = t[i + 1];
            t_c = t[i    ];
            t_w = t[i - 1];
            rhs = 0;
            res = (rhs - a_w * t_w - a_c * t_c - a_e * t_e) / a_c;
            tnew[i] = t_c + res;
            ee_sum_local += res * res;
        }
        
        // dealing with subdomains and boundaries
        if (size > 1) {
            if (rank == 0) {
                // left end, boundary
                #pragma acc serial present(t[0:len], tnew[0:len]) copyin(k, h) copy(ee_sum_local)
                {
                double a_e, a_c, a_w, t_e, t_c, t_w, res, rhs;
                a_e = k / (h * h);
                a_c = - 3.0 * k / (h * h);
                t_e = t[2];
                t_c = t[1];
                rhs = - 2.0 * k * t0 / (h * h);
                res = (rhs - a_c * t_c - a_e * t_e) / a_c;
                tnew[1] = t_c + res;
                ee_sum_local += res * res;
                }
                // right end, subdomain interface
                MPI_Waitall(2, &req[0], MPI_STATUS_IGNORE);
                #pragma acc serial present(t[0:len], tnew[0:len]) copyin(k, h, len) copy(ee_sum_local)
                {
                double a_e, a_c, a_w, t_e, t_c, t_w, res, rhs;
                a_e = k / (h * h);
                a_c = - 2.0 * k / (h * h);
                a_w = k / (h * h);
                t_e = t[len - 1];
                t_c = t[len - 2];
                t_w = t[len - 3];
                rhs = 0;
                res = (rhs - a_w * t_w - a_c * t_c - a_e * t_e) / a_c;
                tnew[len - 2] = t_c + res;
                ee_sum_local += res * res;
                }
            } else if (rank == size - 1) {
                // right end, boundary
                #pragma acc serial present(t[0:len], tnew[0:len]) copyin(k, h, len) copy(ee_sum_local)
                {
                double a_e, a_c, a_w, t_e, t_c, t_w, res, rhs;
                a_c = - 3.0 * k / (h * h);
                a_w = k / (h * h);
                t_c = t[len - 2];
                t_w = t[len - 3];
                rhs = - 2.0 * k * t1 / (h * h);
                res = (rhs - a_w * t_w - a_c * t_c) / a_c;
                tnew[len - 2] = t_c + res;
                ee_sum_local += res * res;
                }
                // left end, subdomain interface
                MPI_Waitall(2, &req[2], MPI_STATUS_IGNORE);
                #pragma acc serial present(t[0:len], tnew[0:len]) copyin(k, h) copy(ee_sum_local)
                {
                double a_e, a_c, a_w, t_e, t_c, t_w, res, rhs;
                a_e = k / (h * h);
                a_c = - 2.0 * k / (h * h);
                a_w = k / (h * h);
                t_e = t[2];
                t_c = t[1];
                t_w = t[0];
                rhs = 0;
                res = (rhs - a_w * t_w - a_c * t_c - a_e * t_e) / a_c;
                tnew[1] = t_c + res;
                ee_sum_local += res * res;
                }
            } else {
                // left end, subdomain interface
                MPI_Waitall(2, &req[2], MPI_STATUS_IGNORE);
                #pragma acc serial present(t[0:len], tnew[0:len]) copyin(k, h) copy(ee_sum_local)
                {
                double a_e, a_c, a_w, t_e, t_c, t_w, res, rhs;
                a_e = k / (h * h);
                a_c = - 2.0 * k / (h * h);
                a_w = k / (h * h);
                t_e = t[2];
                t_c = t[1];
                t_w = t[0];
                rhs = 0;
                res = (rhs - a_w * t_w - a_c * t_c - a_e * t_e) / a_c;
                tnew[1] = t_c + res;
                ee_sum_local += res * res;
                }
                // right end, subdomain interface
                MPI_Waitall(2, &req[0], MPI_STATUS_IGNORE);
                #pragma acc serial present(t[0:len], tnew[0:len]) copyin(k, h, len) copy(ee_sum_local)
                {
                double a_e, a_c, a_w, t_e, t_c, t_w, res, rhs;
                a_e = k / (h * h);
                a_c = - 2.0 * k / (h * h);
                a_w = k / (h * h);
                t_e = t[len - 1];
                t_c = t[len - 2];
                t_w = t[len - 3];
                rhs = 0;
                res = (rhs - a_w * t_w - a_c * t_c - a_e * t_e) / a_c;
                tnew[len - 2] = t_c + res;
                ee_sum_local += res * res;
                }
            }
            MPI_Allreduce(&ee_sum_local, &ee_sum_global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        } else {
            #pragma acc serial present(t[0:len], tnew[0:len]) copyin(k, h, len) copy(ee_sum_local)
            {
            double a_e, a_c, a_w, t_e, t_c, t_w, res, rhs;
            // left end, boundary
            a_e = k / (h * h);
            a_c = -3.0 * k / (h * h);
            t_e = t[2];
            t_c = t[1];
            rhs = - 2.0 * k * t0 / (h * h);
            res = (rhs - a_c * t_c - a_e * t_e) / a_c;
            tnew[1] = t_c + res;
            ee_sum_local += res * res;
            // right end, boundary
            a_c = - 3.0 * k / (h * h);
            a_w = k / (h * h);
            t_c = t[len - 2];
            t_w = t[len - 3];
            rhs = - 2.0 * k * t1 / (h * h);
            res = (rhs - a_w * t_w - a_c * t_c) / a_c;
            tnew[len - 2] = t_c + res;
            ee_sum_local += res * res;
            }
            ee_sum_global = ee_sum_local;
        }

        #pragma acc kernels loop independent present(t[0:len], tnew[0:len]) copyin(len)
        for (int i = 1; i <= len - 2; i ++) {
            t[i] = tnew[i];
        }

        e = sqrt(ee_sum_global / N);
        it ++;
        if (rank == 0 && it % 10000 == 0) {
            printf("(%9d %.5e)\n", it, e);
            fflush(stdout);
        }
    } while (e > 1e-6);
    #pragma acc exit data copyout(t[0:len]) delete(tnew[0:len])
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        end = clock();
        printf("\n");
        printf("%d %.5e %lf\n", it, e, double(end - start) / CLOCKS_PER_SEC);
    }

    MPI_Gather(&t[1], len - 2, MPI_DOUBLE, t_global, len - 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        FILE *fo = fopen("1dlap.csv", "w+t");
        if (fo == nullptr) {
            printf("Err : cannot open output file 1dlap.csv\n");
        } else {
            fprintf(fo, "x,t\n");
            double x = 0.5 * h;
            for (int i = 0; i < N; i ++) {
                fprintf(fo, "%10.3e,%10.3e\n", x, t_global[i]);
                x += h;
            }
            fclose(fo);
        }
    }

    MPI_Finalize();
    return 0;
}
