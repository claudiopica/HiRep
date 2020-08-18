#include <stdlib.h>
#include <stdio.h>
#include "error.h"
#include "logger.h"
#include "data_storage.h"

data_storage_array *allocate_data_storage_array(int n)
{
    data_storage_array *ret = (data_storage_array *)malloc(sizeof(data_storage_array));
    ret->n = n;
    ret->element = (data_storage *)malloc(n * sizeof(data_storage));
    return ret;
}

void allocate_data_storage_element(data_storage_array *cont, int id, int n, int *ni)
{

    error(cont->n <= id, 0, "allocate_data_storage " __FILE__, "wrong id for the storage array element");

    error(n <= 0, 0, "allocate_data_storage " __FILE__, "Number of data dimensions cannot be smaller than zero");

    cont->element[id].n = n;
    cont->element[id].ni = (int *)malloc(n * sizeof(int));

    int size = 1;
    for (int i = 0; i < n; i++)
    {
        error(ni[i] < 1, 0, "allocate_data_storage " __FILE__, "Data sizes cannot be smaller than zero");
        size *= ni[i];
        cont->element[id].ni[i] = ni[i];
    }

    cont->element[id].data = (double *)calloc(size, sizeof(double));
}

void free_data_storage(data_storage_array *dat)
{
    for (int i = 0; i < dat->n; i++)
    {
        free(dat->element[i].data);
        free(dat->element[i].ni);
    }
    free(dat->element);
    free(dat);
}

double *data_storage_element(data_storage_array *dat, int iel, int *ni)
{
    int idx = 0;
    error(iel >= dat->n, 0, "data_storage_element " __FILE__, "Wrong data size");

    for (int i = 0; i < dat->element[iel].n; i++)
    {
        error(ni[i] >= dat->element[iel].ni[i], 0, "data_storage_element " __FILE__, "Wrong data size");
        idx = idx * dat->element[iel].ni[i] + ni[i];
    }
    return dat->element[iel].data + idx;
}

void print_data_storage(data_storage_array *dat)
{
    if (dat == NULL)
    {
        lprintf("WARNING", 0, "Request of printing a not allocated data_storage\n");
        return;
    }
    int size;

    for (int k = 0; k < dat->n; k++)
    {
        lprintf("DATA_STORAGE", 0, "Element %d\n", k);
        size = 1;
        for (int i = 0; i < dat->element[k].n; i++)
            size *= dat->element[k].ni[i];
        int index, i = 0;
        int idx_strg[dat->element[k].n - 1];
        while (i < size)
        {
            index = i / dat->element[k].ni[dat->element[k].n - 1];

            for (int j = dat->element[k].n - 1; j > 0; j--)
            {
                idx_strg[j - 1] = index % dat->element[k].ni[j - 1];
                index /= dat->element[k].ni[j - 1];
            }
            for (int j = 0; j < dat->element[k].n - 1; j++)
                lprintf("DATA_STORAGE", 0, " %d", idx_strg[j]);

            for (int j = 0; j < dat->element[k].ni[dat->element[k].n - 1]; j++, i++)
            {
                lprintf("DATA_STORAGE", 0, " %2.16e", i, size, *(dat->element[k].data + i));
            }
            lprintf("DATA_STORAGE", 0, "\n");
        }
    }
}
