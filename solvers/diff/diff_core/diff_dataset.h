#pragma once

#include "convdiff_core.h"
#include "rom_dataset.h"
#include "hrom_dataset.h"

typedef struct
{
	int    num_ip_each_axis;
	double mat_epsilon;
	int    mat_max_iter;

	double dt;
	double finish_time;
	int    output_interval;

	double* T;
	double* error;
	double* theo_sol;

    /*for ROM input data*/
    double  rom_finish_time;
    int     rom_output_interval;
    int     snapshot_interval;

    int     num_cases;
    double* density_cases;
    double* viscosity_cases;

} VALUES;

typedef struct
{
	const char* directory;

} CONDITIONS;

/*for ROM input data*/
typedef struct
{
    int     num_modes;		    //1領域当たりの最大基底本数
    double  rom_epsilon;	    //基底取得閾値
    int     num_subdomains;     //POD計算領域数
    int     solver_type;
} ROM_PRM;


typedef struct
{
	BBFE_BASIS   basis;
	BBFE_DATA    fe;
	BBFE_BC      bc;
	MONOLIS      monolis;
	MONOLIS_COM  monolis_com;

	CONDITIONS   cond;
	VALUES       vals;

	MONOLIS      monolis0;          // for nonsteady analysis


    /*for ROM*/
	MONOLIS      monolis_rom;
	MONOLIS      monolis_rom0;
	
	VALUES       vals_rom;

	MONOLIS_COM  mono_com_rom;
	MONOLIS_COM  mono_com0;             //初期化関数
	MONOLIS_COM  mono_com_rom_solv;	    //線形ソルバ用コミュニケータ

	ROM_PRM		 rom_prm;		        //input データ用
	ROM		  	 rom;

    HROM		 hrom;
	MONOLIS      monolis_hr;
	MONOLIS      monolis_hr0;

} FE_SYSTEM;