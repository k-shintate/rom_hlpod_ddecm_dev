#pragma once

#include "convdiff_core.h"
#include "nedelec_core.h"
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

	//double* T;
	double* error;
	double* theo_sol;

    double * V;
    double * phi;
    double _Complex * Aphi;
    double* Aphi_time;
    double* Aphi_time_curr;

    double* Aphi_time_rom;
    double* Aphi_time_curr_rom;

    /*for ROM input data*/
    double  rom_finish_time;
    int     rom_output_interval;
    int     snapshot_interval;

    int     num_cases;
    double* density_cases;
    double* viscosity_cases;

	double** v;
	double*  p;

} VALUES;

/*for ROM input data*/
typedef struct
{
    int     num_modes;		    //1領域当たりの最大基底本数
    double  rom_epsilon;	    //基底取得閾値
    int     num_subdomains;     //POD計算領域数
    int     solver_type;
    int     hot_start;
    double  hot_start_time;

} ROM_PRM;

typedef struct
{
	const char* directory;

} CONDITIONS;

typedef struct
{
	BBFE_BASIS   basis;
	BBFE_DATA    fe;
	BBFE_BC      bc;
    BBFE_BC      bc_internal;
    BBFE_BC      bc_NR;
	MONOLIS      monolis;
	MONOLIS      monolis0;
	MONOLIS_COM  monolis_com;
    //MONOLIS_COM  mono_com;

    MONOLIS      monolis_mass;

	CONDITIONS   cond;
	VALUES       vals;

    NEDELEC     ned;

    /*for ROM*/
	MONOLIS      monolis_rom;
	MONOLIS      monolis_rom0;
	
	VALUES       vals_rom;

	MONOLIS_COM  mono_com_rom;
	MONOLIS_COM  mono_com0;	        //初期化関数
	MONOLIS_COM  mono_com_rom_solv;	//線形ソルバ用コミュニケータ

	ROM_PRM		 rom_prm_p;		    //p の input データ用
	ROM_PRM		 rom_prm_v;		    //v の input データ用

	ROM		  	 rom_p;		        //p の snapshot matrix 用関数
	ROM		  	 rom_v;		        //v の snapshot matrix 用関数
	ROM		  	 rom_sups;	        //p+v　の行列演算関連

    VALUES       vals_hrom;

	HROM		 hrom_p;		    //p の snapshot matrix 用関数
	HROM		 hrom_v;		    //v の snapshot matrix 用関数
	HROM		 hrom_sups;	        //p+v　の行列演算関連

	MONOLIS      monolis_hr;
	MONOLIS      monolis_hr0;

} FE_SYSTEM;
