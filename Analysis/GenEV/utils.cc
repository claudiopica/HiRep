
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <cmath>
#include <iomanip> // std::setprecision
#include <cstdlib>
#include <vector>
#include <sstream>
#include <complex>

#include "utils.h"
#include "matrix_algebra.h"
//#include "listrange.h"

static void check_file_data(std::ifstream &stream_file, int *, int *, int *);

static bool data_fine = false;

void check_list_dir_and_data_file(par *apar)
{
    datafile addfile;
    if (apar->Dflag == 1)
    {
        std::ifstream file_with_dirs(apar->list_dir, std::ios::in);
        std::string line;
        if (file_with_dirs.good())
            while (getline(file_with_dirs, line))
            {
                if (!line.empty())
                {
                    if (line.back() != '/')
                        line += "/";

                    addfile.filename = line + apar->vev_name;
                    apar->vevfiles.push_back(addfile);
                }
            }
        else
        {
            std::cerr << "[ERROR][CHECK_LIST_DIR_AND_DATA_FILE] The file \"" << apar->list_dir << "\" does not exist!" << std::endl;
            exit(1);
        }

        file_with_dirs.close();
    }
    else
    {

        addfile.filename = apar->vev_name;
        apar->vevfiles.push_back(addfile);
    }
    std::ifstream data_file;

    apar->numbinjk = 0;

    for (std::vector<datafile>::iterator it_file = apar->vevfiles.begin(); it_file != apar->vevfiles.end(); ++it_file)
    {
        data_file.open(it_file->filename, std::ios::in);
        std::cout << "[INFO][CHECK_LIST_DIR_AND_DATA_FILE] Checking data file \"" << it_file->filename << "\"" << std::endl;
        if (!data_file.good())
        {
            std::cerr << "[ERROR][CHECK_LIST_DIR_AND_DATA_FILE_AND_DATA_FILE] The file \"" << it_file->filename << "\" does not exist!" << std::endl;
            exit(1);
        }
        check_file_data(data_file, &(apar->numop), &(apar->nt), &(it_file->nmeas));

        apar->numbinjk += it_file->nmeas / apar->binwidth;
        data_file.close();
    }

    for (int i = 0; i < apar->numop; i++)
        apar->activeop.push_back(i);

    data_fine = true;
}

void print_parameters(const par *apar)
{
    std::cout << "[INFO][PRINT_PARAMETERS]" << std::endl;
    std::cout << "number of requested corr dist (ndt) ->  " << apar->ndt << std::endl;
    std::cout << "                      vev file name ->  " << apar->vev_name << std::endl;
    std::cout << "                number of operators ->  " << apar->numop << std::endl;
    std::cout << "        number of T valued vev (nt) ->  " << apar->nt << std::endl;
    std::cout << "                   size of each bin ->  " << apar->binwidth << std::endl;
    std::cout << "                     number of bins ->  " << apar->numbinjk << std::endl;
    std::cout << "           inversion point (dT_inv) ->  " << apar->p_inv << std::endl;
    std::cout << "    diagonalization point (dT_diag) ->  " << apar->p_diag << std::endl;
    std::cout << "      Largest allocation size in Gb ->  " << double(apar->numop * apar->numop * apar->ndt * apar->numbinjk * sizeof(double)) / 1073741824 << std::endl;
}

static void check_file_data(std::ifstream &stream_file, int *opcount, int *nidx2, int *nidx1)
{
    std::string line;
    double dword;
    int iword, istart;
    int op_count = 0;
    int nrep = 0;
    int ndata = 1;

    stream_file >> istart;
    stream_file >> iword;

    getline(stream_file, line);
    std::stringstream s(line);

    while (s >> dword)
        op_count++;

    op_count /= 2;

    if (*opcount == 0)
    {
        *opcount = op_count;
    }

    if (*opcount != op_count)
    {
        std::cerr << "[ERROR][check_file_data] Different number of operators " << *opcount << " " << op_count << " between files!" << std::endl;
        exit(1);
    }

    stream_file.seekg(0);

    while (!stream_file.eof())
    {
        stream_file >> iword;

        if (iword != istart)
            break;
        nrep++;
        getline(stream_file, line);
    }

    stream_file.seekg(-sizeof(int), stream_file.cur);
    stream_file >> iword;
    do
    {
        for (int j = 0; j < nrep; j++)
        {
            stream_file >> iword;
            getline(stream_file, line);
            stream_file >> iword;
        }
        ndata++;

        if (stream_file.eof())
            break;
    } while (true);

    if (*nidx2 == 0)
        *nidx2 = nrep;

    if (*nidx2 != nrep)
    {

        std::cerr << "[ERROR][check_file_data] Different number nrep " << *nidx2 << " " << nrep << "(second index) between files!" << std::endl;
        exit(1);
    }

    *nidx1 = ndata;
}

void read_2pt_def(par *apar)
{
    if (!data_fine)
    {
        std::cerr << "Run first check_list_dir_and_data_file to validate the data files " << std::endl;
        exit(1);
    }
    std::string line;
    std::ifstream file_with_defs(apar->cor_def_filename, std::ios::in);
    int p1 = 0, p2 = 0, dt;
    if (file_with_defs.good())
    {
        while (getline(file_with_defs, line))
        {
            if (!line.empty())
            {

                size_t current, previous = 0;
                current = line.find(",");
                while (true)
                {
                    std::string subline = line.substr(previous, current - previous);

                    size_t current2, previous2 = 0;
                    current2 = subline.find("|");
                    dt = 0;
                    while (true)
                    {
                        std::string subsubline = subline.substr(previous2, current2 - previous2);
                        if (sscanf(subsubline.c_str(), "%d-%d", &p1, &p2) == 2)
                        {
                            if (p2 <= p1)
                            {
                                std::cerr << "[ERROR][read_2pt_def] Wrong format in the corellator definition string. \n The correlators are inputed according to the structure x1-y1|x2-y2,x3-y3 with yi>xi and pair separated by the | must have the same distance" << std::endl;
                                exit(1);
                            }
                            if (dt == 0 || dt == p2 - p1)
                            {
                                dt = p2 - p1;
                                apar->corrdef[dt].points.push_back(cpoints(p1, p2));
                            }
                            else
                            {
                                std::cerr << "[ERROR][read_2pt_def] Wrong format in the corellator definition string. \n The correlators are inputed according to the structure x1-y1|x2-y2,x3-y3 with yi>xi and pair separated by the | must have the same distance" << std::endl;
                                exit(1);
                            }
                        }
                        else
                        {
                            std::cerr << "[ERROR][read_2pt_def] Wrong format in the corellator definition string. \n The correlators are inputed according to the structure x1-y1|x2-y2,x3-y3 with yi>xi and pair separated by the | must have the same distance" << std::endl;
                            exit(1);
                        }
                        if (current2 == std::string::npos)
                            break;
                        previous2 = current2 + 1;
                        current2 = subline.find("|", previous2);
                    }

                    if (current == std::string::npos)
                        break;
                    previous = current + 1;
                    current = line.find(",", previous);
                }
            }
        }
        int mindt = 100;

        int idt = 0;
        for (std::map<int, dtcor>::iterator it_dt = apar->corrdef.begin(); it_dt != apar->corrdef.end(); ++it_dt)
        {
            if (it_dt->first < mindt)
                mindt = it_dt->first;
            it_dt->second.index = idt;
            idt++;
        }
        if (apar->xflag == 1)
        {
            if (apar->corrdef.find(apar->p_inv) == apar->corrdef.end())
            {
                std::cerr << "[ERROR][READ_2PT_DEF] The inversion distance is not measured in  " << apar->cor_def_filename << std::endl;
                exit(1);
            }
        }
        else
            apar->p_inv = mindt;

        if (apar->dflag == 1)
        {
            if (apar->corrdef.find(apar->p_diag) == apar->corrdef.end())
            {
                std::cerr << "[ERROR][READ_2PT_DEF] The Diagionalization distance is not measured in  " << apar->cor_def_filename << std::endl;
                exit(1);
            }
        }
        else
            apar->p_diag = apar->p_inv + 1;
        apar->ndt = apar->corrdef.size();
    }
    else
    {
        std::cerr << "[ERROR][READ_2PT_DEF] The file \"" << apar->cor_def_filename << "\" does not exist!" << std::endl;
        exit(1);
    }
}
void read_1pt_map(par *apar)
{
    if (!data_fine)
    {
        std::cerr << "Run first check_list_dir_and_data_file to validate the data files " << std::endl;
        exit(1);
    }
    int istart, iword, nrep = 0;

    std::ifstream data_file;
    std::string line;
    data_file.open(apar->vevfiles[0].filename, std::ios::in);

    data_file >> istart;

    data_file.seekg(0);

    while (!data_file.eof())
    {
        data_file >> iword;
        if (iword != istart)
            break;
        data_file >> iword;
        apar->t2idx[iword] = nrep;
        nrep++;
        getline(data_file, line);
    }

    data_file.close();

    if (apar->t2idx.size() - apar->nt != 0)
    {
        std::cerr << "number of t vev positions differes between read_1pt_map and check_file_data " << std::endl;
        exit(1);
    }
}

void read_data(double *cor_dataout, std::complex<double> *vev_dataout, par *apar)
{

    // dataout -> is the matrix containing the full statistical set: it is created by binning the measurements taken from the files

    int iword;
    double dwordre, dwordim;

    double norm_vev = 1.0 / (apar->binwidth);
    double norm_cor = 1.0 / (apar->binwidth);

    std::ifstream in_file;

    int measurecount = 0;

    std::complex<double> localdata[apar->numop * apar->nt];

    int it1, it2, countdt;
    for (int i = 0; i < apar->numop * apar->numop * apar->ndt * apar->numbinjk; i++)
        cor_dataout[i] = 0.0;

    for (int i = 0; i < apar->numop * apar->nt * apar->numbinjk; i++)
        vev_dataout[i] = 0.0;

    for (std::vector<datafile>::iterator it_file = apar->vevfiles.begin(); it_file != apar->vevfiles.end(); ++it_file)
    {
        std::cout << "[INFO][READ_DATA] Reading data file \"" << it_file->filename << "\"" << std::endl;
        in_file.open(it_file->filename, std::ios::in);

        for (int i = 0; i < it_file->nmeas / apar->binwidth; i++)
        {
            for (int l = 0; l < apar->binwidth; l++)
            {
                for (int j = 0; j < apar->nt; j++)
                {
                    in_file >> iword;
                    in_file >> iword;
                    for (int op1 = 0; op1 < apar->numop; op1++)
                    {
                        in_file >> dwordre;
                        in_file >> dwordim;
                        localdata[op1 + apar->numop * j] = std::complex<double>(dwordre, dwordim);
                        vev_dataout[op1 + apar->numop * (j + apar->nt * measurecount)] += localdata[op1 + apar->numop * j] * norm_vev;
                    }
                }

                countdt = 0;
                for (std::map<int, dtcor>::iterator if_dist = apar->corrdef.begin(); if_dist != apar->corrdef.end(); ++if_dist)
                {
                    norm_cor = 1.0 / (apar->binwidth * if_dist->second.points.size());
                    for (std::vector<cpoints>::iterator pair = if_dist->second.points.begin(); pair != if_dist->second.points.end(); ++pair)
                    {
                        it1 = apar->t2idx[pair->p1];
                        it2 = apar->t2idx[pair->p2];
                        for (int op1 = 0; op1 < apar->numop; op1++)
                            for (int op2 = 0; op2 < apar->numop; op2++)
                            {
                                cor_dataout[op2 + apar->numop * ((op1) + apar->numop * (countdt + measurecount * apar->ndt))] +=
                                    norm_cor * real(conj(localdata[op1 + apar->numop * it1]) * localdata[op2 + apar->numop * it2]);
                            }
                    }
                    countdt++;
                }
            }
            measurecount++;
        }
        in_file.close();
    }
}

void create_simmetrize_corr_jk(par *apar, double *cor_jk, double *cor, std::complex<double> *vev)
{
    double average = 0.0;
    double *vev_avg = new double[apar->numop * apar->nt];
    double *vev_jk = new double[apar->numop * apar->nt * apar->numbinjk];
    if (apar->vflag == 1)
    {

        std::cout << "[INFO][SIMMETRIZE_CORR_JK] Subtracting the vev from two point functions\n";
        for (int op1 = 0; op1 < apar->numop; op1++)
            for (int nt1 = 0; nt1 < apar->nt; nt1++)
            {
                average = 0.0;

                for (int bin = 0; bin < apar->numbinjk; bin++)
                    average += vev[op1 + apar->numop * (nt1 + apar->nt * bin)].real();

                vev_avg[op1 + apar->numop * nt1] = average / apar->numbinjk;

                for (int bin = 0; bin < apar->numbinjk; bin++)
                    vev_jk[op1 + apar->numop * (nt1 + apar->nt * bin)] =
                        (average - vev[op1 + apar->numop * (nt1 + apar->nt * bin)].real()) / (apar->numbinjk - 1);
            }
    }

    int dtcounter = 0;
    std::cout << "[INFO][SIMMETRIZE_CORR_JK] Creating the JK bins\n";
    int it1, it2;
    for (int op1 = 0; op1 < apar->numop; op1++)
        for (int op2 = 0; op2 < apar->numop; op2++)
        {
            dtcounter = 0;

            for (std::map<int, dtcor>::iterator if_dist = apar->corrdef.begin(); if_dist != apar->corrdef.end(); ++if_dist)
            {

                average = 0.0;
                for (int bin = 0; bin < apar->numbinjk; bin++)
                {
                    average += cor[op2 + apar->numop * ((op1) + apar->numop * (dtcounter + bin * apar->ndt))];
                }

                for (int bin = 0; bin < apar->numbinjk; bin++)
                {
                    cor_jk[op2 + apar->numop * ((op1) + apar->numop * (dtcounter + bin * apar->ndt))] =
                        (average - cor[op2 + apar->numop * ((op1) + apar->numop * (dtcounter + bin * apar->ndt))]) / (apar->numbinjk - 1.0);
                    if (apar->vflag == 1)
                    {
                        for (std::vector<cpoints>::iterator pair = if_dist->second.points.begin(); pair != if_dist->second.points.end(); ++pair)
                        {
                            it1 = apar->t2idx[pair->p1];
                            it2 = apar->t2idx[pair->p2];
                            cor_jk[op2 + apar->numop * ((op1) + apar->numop * (dtcounter + bin * apar->ndt))] -=
                                vev_jk[op2 + apar->numop * (it2 + apar->nt * bin)] * vev_jk[op1 + apar->numop * (it1 + apar->nt * bin)] / if_dist->second.points.size();
                        }
                    }
                }

                dtcounter++;
            }
        }

    std::cout << "[INFO][SIMMETRIZE_CORR_JK] Simmetrizing the correlation matrix" << std::endl;

    double norm = 0.5;

    for (int op1 = 0; op1 < apar->numop; op1++)
    {

        for (int op2 = op1; op2 < apar->numop; op2++)
        {
            for (int t = 0; t < apar->ndt; t++)
            {
                for (int bin = 0; bin < apar->numbinjk; bin++)
                {
                    cor_jk[op2 + apar->numop * ((op1) + apar->numop * (t + bin * apar->ndt))] = norm *
                                                                                                (cor_jk[op2 + apar->numop * ((op1) + apar->numop * (t + bin * apar->ndt))] +
                                                                                                 cor_jk[op1 + apar->numop * ((op2) + apar->numop * (t + bin * apar->ndt))]);
                    cor_jk[op1 + apar->numop * ((op2) + apar->numop * (t + bin * apar->ndt))] = cor_jk[op2 + apar->numop * ((op1) + apar->numop * (t + bin * apar->ndt))];
                }
            }
        }
    }
    delete[] vev_avg;
    delete[] vev_jk;
}

void deactivate_pull_corr_jk(par *apar, double *cor_jk, double pull, int dt)
{

    std::cout << "[INFO][deactivate_pull_corr_jk] Deactivating the op with pull<" << pull << " at dt=" << dt << std::endl;
    int idt = apar->corrdef[dt].index;

    for (int op1 = 0; op1 < apar->numop; op1++)
    {
        double corjk_var = 0., corjk_avg = 0.0;

        for (int bin = 0; bin < apar->numbinjk; bin++)
        {
            corjk_avg += cor_jk[op1 + apar->numop * ((op1) + apar->numop * (idt + bin * apar->ndt))];
            corjk_var += cor_jk[op1 + apar->numop * ((op1) + apar->numop * (idt + bin * apar->ndt))] *
                         cor_jk[op1 + apar->numop * ((op1) + apar->numop * (idt + bin * apar->ndt))];
        }
        corjk_avg /= apar->numbinjk;
        corjk_var /= apar->numbinjk;
        corjk_var = sqrt((double)(apar->numbinjk - 1.0) * (corjk_var - corjk_avg * corjk_avg));
        if (corjk_avg / corjk_var < pull)
            deactivate_op(apar, op1);
    }
}

void print_diag_cor_avg(par *apar, int nop, double *cor_jk)
{

    std::cout << "[INFO][print_diag_cor_avg] Diagonal correlators jk analysis" << std::endl;
    int t;
    double corjk_var = 0., corjk_avg = 0.0;

    for (std::map<int, dtcor>::iterator if_dist = apar->corrdef.begin(); if_dist != apar->corrdef.end(); ++if_dist)
    {
        t = if_dist->second.index;
        std::cout << if_dist->first << " ";

        for (int op1 = 0; op1 < (int)(nop); op1++)
        {

            corjk_var = 0., corjk_avg = 0.0;

            for (int bin = 0; bin < apar->numbinjk; bin++)
            {
                corjk_avg += cor_jk[op1 + nop * ((op1) + nop * (t + bin * apar->ndt))];
                corjk_var += cor_jk[op1 + nop * ((op1) + nop * (t + bin * apar->ndt))] *
                             cor_jk[op1 + nop * ((op1) + nop * (t + bin * apar->ndt))];
            }
            corjk_avg /= apar->numbinjk;
            corjk_var /= apar->numbinjk;
            corjk_var = sqrt((double)(apar->numbinjk - 1.0) * (corjk_var - corjk_avg * corjk_avg));
            std::cout << corjk_avg << " " << corjk_var << " ";
        }
        std::cout << std::endl;
    }
}

void deactivate_op(par *apar, int opid)
{
    if (opid > apar->numop)
    {
        std::cout << "[WARNING][deactivate_op] The requested operator index is larger than the available list" << std::endl;
        return;
    }

    for (int i = 0; i < (int)(apar->activeop.size()); i++)
        if (opid == apar->activeop[i])
        {
            std::cout << "[INFO][deactivate_op] Deactivating op=" << opid << std::endl;
            apar->activeop.erase(apar->activeop.begin() + i);
            return;
        }
    std::cout << "[WARNING][deactivate_op] The operator " << opid << " is already not active" << std::endl;
}

void normalize_corr_jk(par *apar, double *cor_jk, int dt)
{
    std::cout << "[INFO][NORMALIZE_CORR_JK] Normalizing at dt=" << apar->corrdef.begin()->first << std::endl;
    int idt = apar->corrdef[dt].index;

    double *corjk_avg = new double[apar->numop];

    for (int op1 = 0; op1 < apar->numop; op1++)
    {
        corjk_avg[op1] = 0.;

        for (int bin = 0; bin < apar->numbinjk; bin++)
        {
            corjk_avg[op1] += cor_jk[op1 + apar->numop * ((op1) + apar->numop * (idt + bin * apar->ndt))];
        }
        corjk_avg[op1] /= apar->numbinjk;
    }
    double norm;

    for (int op1 = 0; op1 < apar->numop; op1++)
        for (int op2 = 0; op2 < apar->numop; op2++)
        {
            norm = 1.0 / sqrt(abs(corjk_avg[op1] * corjk_avg[op2]));
            for (int t = 0; t < apar->ndt; t++)
            {
                for (int bin = 0; bin < apar->numbinjk; bin++)
                {
                    cor_jk[op1 + apar->numop * ((op2) + apar->numop * (t + bin * apar->ndt))] *= norm;
                }
            }
        }
}

double *generate_cut_matrix(par *apar, double *cor_jk)
{
    double *retmat = new double[apar->activeop.size() * apar->activeop.size() * apar->ndt * apar->numbinjk];
    int op1 = 0, op2 = 0;
    for (std::vector<int>::iterator it_op1 = apar->activeop.begin(); it_op1 != apar->activeop.end(); ++it_op1, ++op1)
    {
        op2 = 0;
        for (std::vector<int>::iterator it_op2 = apar->activeop.begin(); it_op2 != apar->activeop.end(); ++it_op2, ++op2)
        {
            for (int idt = 0; idt < apar->ndt; idt++)
                for (int bin = 0; bin < apar->numbinjk; bin++)
                {
                    retmat[op1 + apar->activeop.size() * ((op2) + apar->activeop.size() * (idt + bin * apar->ndt))] =
                        cor_jk[*it_op1 + apar->numop * ((*it_op2) + apar->numop * (idt + bin * apar->ndt))];
                }
        }
    }
    return retmat;
}

void get_cm1c_matrix(const par *apar, double *jk_mat, int jkstep, double *mat_inverse, int point, double *jk_corr_mat_t)
{
    int mat_size = apar->activeop.size();

    double *mat = jk_mat + mat_size * mat_size * (point + (jkstep)*apar->ndt);

    for (int j = 0; j < mat_size; j++)
        for (int i = 0; i < mat_size; i++)
        {
            jk_corr_mat_t[i * mat_size + j] = 0.0;
            for (int k = 0; k < mat_size; k++)
                jk_corr_mat_t[i * mat_size + j] += mat_inverse[i * mat_size + k] * mat[k * mat_size + j];
        }
}