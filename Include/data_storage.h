
#ifndef DATA_STORAGE_H
#define DATA_STORAGE_H
typedef enum
{
    STORE,
    DONTSTORE
} storage_switch;

typedef struct _data_storage
{
    int n;
    int *ni;
    double *data;
} data_storage;

typedef struct _data_storage_array
{
    int n;
    data_storage *element;

} data_storage_array;

data_storage_array *allocate_data_storage_array(int n);

void allocate_data_storage_element(data_storage_array *cont, int id, int n, int *ni);

void free_data_storage(data_storage_array *dat);

double *data_storage_element(data_storage_array *dat, int iel, int *ni);

void print_data_storage(data_storage_array *dat);
#endif