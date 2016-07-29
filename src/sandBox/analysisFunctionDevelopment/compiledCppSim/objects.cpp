
#include "objects.h"
#include "synapses_classes.h"
#include "brianlib/clocks.h"
#include "brianlib/dynamic_array.h"
#include "brianlib/stdint_compat.h"
#include "network.h"
#include "randomkit.h"
#include<vector>
#include<iostream>
#include<fstream>

namespace brian {

std::vector< rk_state* > _mersenne_twister_states;

//////////////// networks /////////////////
Network magicnetwork;

//////////////// arrays ///////////////////
double * _array_defaultclock_dt;
const int _num__array_defaultclock_dt = 1;
double * _array_defaultclock_t;
const int _num__array_defaultclock_t = 1;
uint64_t * _array_defaultclock_timestep;
const int _num__array_defaultclock_timestep = 1;
int32_t * _array_neurongroup__spikespace;
const int _num__array_neurongroup__spikespace = 2;
int32_t * _array_neurongroup_i;
const int _num__array_neurongroup_i = 1;
double * _array_neurongroup_I;
const int _num__array_neurongroup_I = 1;
double * _array_neurongroup_L;
const int _num__array_neurongroup_L = 1;
double * _array_neurongroup_lastspike;
const int _num__array_neurongroup_lastspike = 1;
bool * _array_neurongroup_not_refractory;
const int _num__array_neurongroup_not_refractory = 1;
double * _array_neurongroup_u;
const int _num__array_neurongroup_u = 1;
double * _array_neurongroup_v;
const int _num__array_neurongroup_v = 1;
int32_t * _array_poissongroup__spikespace;
const int _num__array_poissongroup__spikespace = 101;
int32_t * _array_poissongroup_i;
const int _num__array_poissongroup_i = 100;
double * _array_poissongroup_rates;
const int _num__array_poissongroup_rates = 100;
int32_t * _array_spikemonitor_1__source_idx;
const int _num__array_spikemonitor_1__source_idx = 1;
int32_t * _array_spikemonitor_1_count;
const int _num__array_spikemonitor_1_count = 1;
int32_t * _array_spikemonitor_1_N;
const int _num__array_spikemonitor_1_N = 1;
int32_t * _array_spikemonitor__source_idx;
const int _num__array_spikemonitor__source_idx = 100;
int32_t * _array_spikemonitor_count;
const int _num__array_spikemonitor_count = 100;
int32_t * _array_spikemonitor_N;
const int _num__array_spikemonitor_N = 1;
int32_t * _array_statemonitor_1__indices;
const int _num__array_statemonitor_1__indices = 1;
double * _array_statemonitor_1_clock_dt;
const int _num__array_statemonitor_1_clock_dt = 1;
double * _array_statemonitor_1_clock_t;
const int _num__array_statemonitor_1_clock_t = 1;
uint64_t * _array_statemonitor_1_clock_timestep;
const int _num__array_statemonitor_1_clock_timestep = 1;
double * _array_statemonitor_1_L;
const int _num__array_statemonitor_1_L = (0, 1);
int32_t * _array_statemonitor_1_N;
const int _num__array_statemonitor_1_N = 1;
double * _array_statemonitor_1_u;
const int _num__array_statemonitor_1_u = (0, 1);
double * _array_statemonitor_1_v;
const int _num__array_statemonitor_1_v = (0, 1);
int32_t * _array_statemonitor__indices;
const int _num__array_statemonitor__indices = 100;
double * _array_statemonitor_clock_dt;
const int _num__array_statemonitor_clock_dt = 1;
double * _array_statemonitor_clock_t;
const int _num__array_statemonitor_clock_t = 1;
uint64_t * _array_statemonitor_clock_timestep;
const int _num__array_statemonitor_clock_timestep = 1;
double * _array_statemonitor_K;
const int _num__array_statemonitor_K = (0, 100);
int32_t * _array_statemonitor_N;
const int _num__array_statemonitor_N = 1;
double * _array_statemonitor_R_hat;
const int _num__array_statemonitor_R_hat = (0, 100);
double * _array_statemonitor_w;
const int _num__array_statemonitor_w = (0, 100);
int32_t * _array_synapses_N;
const int _num__array_synapses_N = 1;

//////////////// dynamic arrays 1d /////////
std::vector<int32_t> _dynamic_array_spikemonitor_1_i;
std::vector<double> _dynamic_array_spikemonitor_1_t;
std::vector<int32_t> _dynamic_array_spikemonitor_i;
std::vector<double> _dynamic_array_spikemonitor_t;
std::vector<double> _dynamic_array_statemonitor_1_t;
std::vector<double> _dynamic_array_statemonitor_t;
std::vector<int32_t> _dynamic_array_synapses__synaptic_post;
std::vector<int32_t> _dynamic_array_synapses__synaptic_pre;
std::vector<double> _dynamic_array_synapses_apost;
std::vector<double> _dynamic_array_synapses_apre;
std::vector<double> _dynamic_array_synapses_delay;
std::vector<double> _dynamic_array_synapses_delay_1;
std::vector<double> _dynamic_array_synapses_lastupdate;
std::vector<int32_t> _dynamic_array_synapses_N_incoming;
std::vector<int32_t> _dynamic_array_synapses_N_outgoing;
std::vector<double> _dynamic_array_synapses_R_hat;
std::vector<double> _dynamic_array_synapses_w;

//////////////// dynamic arrays 2d /////////
DynamicArray2D<double> _dynamic_array_statemonitor_1_L;
DynamicArray2D<double> _dynamic_array_statemonitor_1_u;
DynamicArray2D<double> _dynamic_array_statemonitor_1_v;
DynamicArray2D<double> _dynamic_array_statemonitor_K;
DynamicArray2D<double> _dynamic_array_statemonitor_R_hat;
DynamicArray2D<double> _dynamic_array_statemonitor_w;

/////////////// static arrays /////////////
double * _static_array__array_poissongroup_rates;
const int _num__static_array__array_poissongroup_rates = 100;
int32_t * _static_array__array_statemonitor__indices;
const int _num__static_array__array_statemonitor__indices = 100;

//////////////// synapses /////////////////
// synapses
SynapticPathway<double> synapses_post(
		_dynamic_array_synapses_delay_1,
		_dynamic_array_synapses__synaptic_post,
		0, 1);
SynapticPathway<double> synapses_pre(
		_dynamic_array_synapses_delay,
		_dynamic_array_synapses__synaptic_pre,
		0, 100);

//////////////// clocks ///////////////////
Clock defaultclock;  // attributes will be set in run.cpp
Clock statemonitor_1_clock;  // attributes will be set in run.cpp
Clock statemonitor_clock;  // attributes will be set in run.cpp

// Profiling information for each code object
double neurongroup_resetter_codeobject_profiling_info = 0.0;
double neurongroup_stateupdater_codeobject_profiling_info = 0.0;
double neurongroup_thresholder_codeobject_profiling_info = 0.0;
double poissongroup_thresholder_codeobject_profiling_info = 0.0;
double spikemonitor_1_codeobject_profiling_info = 0.0;
double spikemonitor_codeobject_profiling_info = 0.0;
double statemonitor_1_codeobject_profiling_info = 0.0;
double statemonitor_codeobject_profiling_info = 0.0;
double synapses_post_codeobject_profiling_info = 0.0;
double synapses_post_initialise_queue_profiling_info = 0.0;
double synapses_post_push_spikes_profiling_info = 0.0;
double synapses_pre_codeobject_profiling_info = 0.0;
double synapses_pre_initialise_queue_profiling_info = 0.0;
double synapses_pre_push_spikes_profiling_info = 0.0;
double synapses_stateupdater_codeobject_profiling_info = 0.0;
double synapses_synapses_create_generator_codeobject_profiling_info = 0.0;

}

void _init_arrays()
{
	using namespace brian;

    // Arrays initialized to 0
	_array_statemonitor__indices = new int32_t[100];
	
	for(int i=0; i<100; i++) _array_statemonitor__indices[i] = 0;
	_array_statemonitor_1__indices = new int32_t[1];
	
	for(int i=0; i<1; i++) _array_statemonitor_1__indices[i] = 0;
	_array_spikemonitor__source_idx = new int32_t[100];
	
	for(int i=0; i<100; i++) _array_spikemonitor__source_idx[i] = 0;
	_array_spikemonitor_1__source_idx = new int32_t[1];
	
	for(int i=0; i<1; i++) _array_spikemonitor_1__source_idx[i] = 0;
	_array_neurongroup__spikespace = new int32_t[2];
	
	for(int i=0; i<2; i++) _array_neurongroup__spikespace[i] = 0;
	_array_poissongroup__spikespace = new int32_t[101];
	
	for(int i=0; i<101; i++) _array_poissongroup__spikespace[i] = 0;
	_array_spikemonitor_count = new int32_t[100];
	
	for(int i=0; i<100; i++) _array_spikemonitor_count[i] = 0;
	_array_spikemonitor_1_count = new int32_t[1];
	
	for(int i=0; i<1; i++) _array_spikemonitor_1_count[i] = 0;
	_array_defaultclock_dt = new double[1];
	
	for(int i=0; i<1; i++) _array_defaultclock_dt[i] = 0;
	_array_statemonitor_clock_dt = new double[1];
	
	for(int i=0; i<1; i++) _array_statemonitor_clock_dt[i] = 0;
	_array_statemonitor_1_clock_dt = new double[1];
	
	for(int i=0; i<1; i++) _array_statemonitor_1_clock_dt[i] = 0;
	_array_neurongroup_i = new int32_t[1];
	
	for(int i=0; i<1; i++) _array_neurongroup_i[i] = 0;
	_array_neurongroup_I = new double[1];
	
	for(int i=0; i<1; i++) _array_neurongroup_I[i] = 0;
	_array_poissongroup_i = new int32_t[100];
	
	for(int i=0; i<100; i++) _array_poissongroup_i[i] = 0;
	_array_neurongroup_L = new double[1];
	
	for(int i=0; i<1; i++) _array_neurongroup_L[i] = 0;
	_array_neurongroup_lastspike = new double[1];
	
	for(int i=0; i<1; i++) _array_neurongroup_lastspike[i] = 0;
	_array_synapses_N = new int32_t[1];
	
	for(int i=0; i<1; i++) _array_synapses_N[i] = 0;
	_array_statemonitor_N = new int32_t[1];
	
	for(int i=0; i<1; i++) _array_statemonitor_N[i] = 0;
	_array_spikemonitor_N = new int32_t[1];
	
	for(int i=0; i<1; i++) _array_spikemonitor_N[i] = 0;
	_array_spikemonitor_1_N = new int32_t[1];
	
	for(int i=0; i<1; i++) _array_spikemonitor_1_N[i] = 0;
	_array_statemonitor_1_N = new int32_t[1];
	
	for(int i=0; i<1; i++) _array_statemonitor_1_N[i] = 0;
	_array_neurongroup_not_refractory = new bool[1];
	
	for(int i=0; i<1; i++) _array_neurongroup_not_refractory[i] = 0;
	_array_poissongroup_rates = new double[100];
	
	for(int i=0; i<100; i++) _array_poissongroup_rates[i] = 0;
	_array_defaultclock_t = new double[1];
	
	for(int i=0; i<1; i++) _array_defaultclock_t[i] = 0;
	_array_statemonitor_clock_t = new double[1];
	
	for(int i=0; i<1; i++) _array_statemonitor_clock_t[i] = 0;
	_array_statemonitor_1_clock_t = new double[1];
	
	for(int i=0; i<1; i++) _array_statemonitor_1_clock_t[i] = 0;
	_array_defaultclock_timestep = new uint64_t[1];
	
	for(int i=0; i<1; i++) _array_defaultclock_timestep[i] = 0;
	_array_statemonitor_clock_timestep = new uint64_t[1];
	
	for(int i=0; i<1; i++) _array_statemonitor_clock_timestep[i] = 0;
	_array_statemonitor_1_clock_timestep = new uint64_t[1];
	
	for(int i=0; i<1; i++) _array_statemonitor_1_clock_timestep[i] = 0;
	_array_neurongroup_u = new double[1];
	
	for(int i=0; i<1; i++) _array_neurongroup_u[i] = 0;
	_array_neurongroup_v = new double[1];
	
	for(int i=0; i<1; i++) _array_neurongroup_v[i] = 0;

	// Arrays initialized to an "arange"
	_array_spikemonitor__source_idx = new int32_t[100];
	
	for(int i=0; i<100; i++) _array_spikemonitor__source_idx[i] = 0 + i;
	_array_spikemonitor_1__source_idx = new int32_t[1];
	
	for(int i=0; i<1; i++) _array_spikemonitor_1__source_idx[i] = 0 + i;
	_array_neurongroup_i = new int32_t[1];
	
	for(int i=0; i<1; i++) _array_neurongroup_i[i] = 0 + i;
	_array_poissongroup_i = new int32_t[100];
	
	for(int i=0; i<100; i++) _array_poissongroup_i[i] = 0 + i;

	// static arrays
	_static_array__array_poissongroup_rates = new double[100];
	_static_array__array_statemonitor__indices = new int32_t[100];

	// Random number generator states
	for (int i=0; i<1; i++)
	    _mersenne_twister_states.push_back(new rk_state());
}

void _load_arrays()
{
	using namespace brian;

	ifstream f_static_array__array_poissongroup_rates;
	f_static_array__array_poissongroup_rates.open("static_arrays/_static_array__array_poissongroup_rates", ios::in | ios::binary);
	if(f_static_array__array_poissongroup_rates.is_open())
	{
		f_static_array__array_poissongroup_rates.read(reinterpret_cast<char*>(_static_array__array_poissongroup_rates), 100*sizeof(double));
	} else
	{
		std::cout << "Error opening static array _static_array__array_poissongroup_rates." << endl;
	}
	ifstream f_static_array__array_statemonitor__indices;
	f_static_array__array_statemonitor__indices.open("static_arrays/_static_array__array_statemonitor__indices", ios::in | ios::binary);
	if(f_static_array__array_statemonitor__indices.is_open())
	{
		f_static_array__array_statemonitor__indices.read(reinterpret_cast<char*>(_static_array__array_statemonitor__indices), 100*sizeof(int32_t));
	} else
	{
		std::cout << "Error opening static array _static_array__array_statemonitor__indices." << endl;
	}
}

void _write_arrays()
{
	using namespace brian;

	ofstream outfile__array_defaultclock_dt;
	outfile__array_defaultclock_dt.open("results\\_array_defaultclock_dt_614249772712282910", ios::binary | ios::out);
	if(outfile__array_defaultclock_dt.is_open())
	{
		outfile__array_defaultclock_dt.write(reinterpret_cast<char*>(_array_defaultclock_dt), 1*sizeof(_array_defaultclock_dt[0]));
		outfile__array_defaultclock_dt.close();
	} else
	{
		std::cout << "Error writing output file for _array_defaultclock_dt." << endl;
	}
	ofstream outfile__array_defaultclock_t;
	outfile__array_defaultclock_t.open("results\\_array_defaultclock_t_-30567324951851406", ios::binary | ios::out);
	if(outfile__array_defaultclock_t.is_open())
	{
		outfile__array_defaultclock_t.write(reinterpret_cast<char*>(_array_defaultclock_t), 1*sizeof(_array_defaultclock_t[0]));
		outfile__array_defaultclock_t.close();
	} else
	{
		std::cout << "Error writing output file for _array_defaultclock_t." << endl;
	}
	ofstream outfile__array_defaultclock_timestep;
	outfile__array_defaultclock_timestep.open("results\\_array_defaultclock_timestep_6356348833222527003", ios::binary | ios::out);
	if(outfile__array_defaultclock_timestep.is_open())
	{
		outfile__array_defaultclock_timestep.write(reinterpret_cast<char*>(_array_defaultclock_timestep), 1*sizeof(_array_defaultclock_timestep[0]));
		outfile__array_defaultclock_timestep.close();
	} else
	{
		std::cout << "Error writing output file for _array_defaultclock_timestep." << endl;
	}
	ofstream outfile__array_neurongroup__spikespace;
	outfile__array_neurongroup__spikespace.open("results\\_array_neurongroup__spikespace_-2370811854520053979", ios::binary | ios::out);
	if(outfile__array_neurongroup__spikespace.is_open())
	{
		outfile__array_neurongroup__spikespace.write(reinterpret_cast<char*>(_array_neurongroup__spikespace), 2*sizeof(_array_neurongroup__spikespace[0]));
		outfile__array_neurongroup__spikespace.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup__spikespace." << endl;
	}
	ofstream outfile__array_neurongroup_i;
	outfile__array_neurongroup_i.open("results\\_array_neurongroup_i_5340428886587075641", ios::binary | ios::out);
	if(outfile__array_neurongroup_i.is_open())
	{
		outfile__array_neurongroup_i.write(reinterpret_cast<char*>(_array_neurongroup_i), 1*sizeof(_array_neurongroup_i[0]));
		outfile__array_neurongroup_i.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_i." << endl;
	}
	ofstream outfile__array_neurongroup_I;
	outfile__array_neurongroup_I.open("results\\_array_neurongroup_I_-5106640032409316959", ios::binary | ios::out);
	if(outfile__array_neurongroup_I.is_open())
	{
		outfile__array_neurongroup_I.write(reinterpret_cast<char*>(_array_neurongroup_I), 1*sizeof(_array_neurongroup_I[0]));
		outfile__array_neurongroup_I.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_I." << endl;
	}
	ofstream outfile__array_neurongroup_L;
	outfile__array_neurongroup_L.open("results\\_array_neurongroup_L_-6509986600802519240", ios::binary | ios::out);
	if(outfile__array_neurongroup_L.is_open())
	{
		outfile__array_neurongroup_L.write(reinterpret_cast<char*>(_array_neurongroup_L), 1*sizeof(_array_neurongroup_L[0]));
		outfile__array_neurongroup_L.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_L." << endl;
	}
	ofstream outfile__array_neurongroup_lastspike;
	outfile__array_neurongroup_lastspike.open("results\\_array_neurongroup_lastspike_-1121403566012987508", ios::binary | ios::out);
	if(outfile__array_neurongroup_lastspike.is_open())
	{
		outfile__array_neurongroup_lastspike.write(reinterpret_cast<char*>(_array_neurongroup_lastspike), 1*sizeof(_array_neurongroup_lastspike[0]));
		outfile__array_neurongroup_lastspike.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_lastspike." << endl;
	}
	ofstream outfile__array_neurongroup_not_refractory;
	outfile__array_neurongroup_not_refractory.open("results\\_array_neurongroup_not_refractory_5597795296270559712", ios::binary | ios::out);
	if(outfile__array_neurongroup_not_refractory.is_open())
	{
		outfile__array_neurongroup_not_refractory.write(reinterpret_cast<char*>(_array_neurongroup_not_refractory), 1*sizeof(_array_neurongroup_not_refractory[0]));
		outfile__array_neurongroup_not_refractory.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_not_refractory." << endl;
	}
	ofstream outfile__array_neurongroup_u;
	outfile__array_neurongroup_u.open("results\\_array_neurongroup_u_-8656701547261747352", ios::binary | ios::out);
	if(outfile__array_neurongroup_u.is_open())
	{
		outfile__array_neurongroup_u.write(reinterpret_cast<char*>(_array_neurongroup_u), 1*sizeof(_array_neurongroup_u[0]));
		outfile__array_neurongroup_u.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_u." << endl;
	}
	ofstream outfile__array_neurongroup_v;
	outfile__array_neurongroup_v.open("results\\_array_neurongroup_v_-1114623716179288880", ios::binary | ios::out);
	if(outfile__array_neurongroup_v.is_open())
	{
		outfile__array_neurongroup_v.write(reinterpret_cast<char*>(_array_neurongroup_v), 1*sizeof(_array_neurongroup_v[0]));
		outfile__array_neurongroup_v.close();
	} else
	{
		std::cout << "Error writing output file for _array_neurongroup_v." << endl;
	}
	ofstream outfile__array_poissongroup__spikespace;
	outfile__array_poissongroup__spikespace.open("results\\_array_poissongroup__spikespace_7109315386570648310", ios::binary | ios::out);
	if(outfile__array_poissongroup__spikespace.is_open())
	{
		outfile__array_poissongroup__spikespace.write(reinterpret_cast<char*>(_array_poissongroup__spikespace), 101*sizeof(_array_poissongroup__spikespace[0]));
		outfile__array_poissongroup__spikespace.close();
	} else
	{
		std::cout << "Error writing output file for _array_poissongroup__spikespace." << endl;
	}
	ofstream outfile__array_poissongroup_i;
	outfile__array_poissongroup_i.open("results\\_array_poissongroup_i_-7372998049886180480", ios::binary | ios::out);
	if(outfile__array_poissongroup_i.is_open())
	{
		outfile__array_poissongroup_i.write(reinterpret_cast<char*>(_array_poissongroup_i), 100*sizeof(_array_poissongroup_i[0]));
		outfile__array_poissongroup_i.close();
	} else
	{
		std::cout << "Error writing output file for _array_poissongroup_i." << endl;
	}
	ofstream outfile__array_poissongroup_rates;
	outfile__array_poissongroup_rates.open("results\\_array_poissongroup_rates_2888448766509021741", ios::binary | ios::out);
	if(outfile__array_poissongroup_rates.is_open())
	{
		outfile__array_poissongroup_rates.write(reinterpret_cast<char*>(_array_poissongroup_rates), 100*sizeof(_array_poissongroup_rates[0]));
		outfile__array_poissongroup_rates.close();
	} else
	{
		std::cout << "Error writing output file for _array_poissongroup_rates." << endl;
	}
	ofstream outfile__array_spikemonitor_1__source_idx;
	outfile__array_spikemonitor_1__source_idx.open("results\\_array_spikemonitor_1__source_idx_-6579324185992451374", ios::binary | ios::out);
	if(outfile__array_spikemonitor_1__source_idx.is_open())
	{
		outfile__array_spikemonitor_1__source_idx.write(reinterpret_cast<char*>(_array_spikemonitor_1__source_idx), 1*sizeof(_array_spikemonitor_1__source_idx[0]));
		outfile__array_spikemonitor_1__source_idx.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor_1__source_idx." << endl;
	}
	ofstream outfile__array_spikemonitor_1_count;
	outfile__array_spikemonitor_1_count.open("results\\_array_spikemonitor_1_count_-4692669552264718276", ios::binary | ios::out);
	if(outfile__array_spikemonitor_1_count.is_open())
	{
		outfile__array_spikemonitor_1_count.write(reinterpret_cast<char*>(_array_spikemonitor_1_count), 1*sizeof(_array_spikemonitor_1_count[0]));
		outfile__array_spikemonitor_1_count.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor_1_count." << endl;
	}
	ofstream outfile__array_spikemonitor_1_N;
	outfile__array_spikemonitor_1_N.open("results\\_array_spikemonitor_1_N_7756002828134391085", ios::binary | ios::out);
	if(outfile__array_spikemonitor_1_N.is_open())
	{
		outfile__array_spikemonitor_1_N.write(reinterpret_cast<char*>(_array_spikemonitor_1_N), 1*sizeof(_array_spikemonitor_1_N[0]));
		outfile__array_spikemonitor_1_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor_1_N." << endl;
	}
	ofstream outfile__array_spikemonitor__source_idx;
	outfile__array_spikemonitor__source_idx.open("results\\_array_spikemonitor__source_idx_4638487192180297482", ios::binary | ios::out);
	if(outfile__array_spikemonitor__source_idx.is_open())
	{
		outfile__array_spikemonitor__source_idx.write(reinterpret_cast<char*>(_array_spikemonitor__source_idx), 100*sizeof(_array_spikemonitor__source_idx[0]));
		outfile__array_spikemonitor__source_idx.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor__source_idx." << endl;
	}
	ofstream outfile__array_spikemonitor_count;
	outfile__array_spikemonitor_count.open("results\\_array_spikemonitor_count_-4416897708574000428", ios::binary | ios::out);
	if(outfile__array_spikemonitor_count.is_open())
	{
		outfile__array_spikemonitor_count.write(reinterpret_cast<char*>(_array_spikemonitor_count), 100*sizeof(_array_spikemonitor_count[0]));
		outfile__array_spikemonitor_count.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor_count." << endl;
	}
	ofstream outfile__array_spikemonitor_N;
	outfile__array_spikemonitor_N.open("results\\_array_spikemonitor_N_7019897520583539084", ios::binary | ios::out);
	if(outfile__array_spikemonitor_N.is_open())
	{
		outfile__array_spikemonitor_N.write(reinterpret_cast<char*>(_array_spikemonitor_N), 1*sizeof(_array_spikemonitor_N[0]));
		outfile__array_spikemonitor_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_spikemonitor_N." << endl;
	}
	ofstream outfile__array_statemonitor_1__indices;
	outfile__array_statemonitor_1__indices.open("results\\_array_statemonitor_1__indices_-7918381090160206349", ios::binary | ios::out);
	if(outfile__array_statemonitor_1__indices.is_open())
	{
		outfile__array_statemonitor_1__indices.write(reinterpret_cast<char*>(_array_statemonitor_1__indices), 1*sizeof(_array_statemonitor_1__indices[0]));
		outfile__array_statemonitor_1__indices.close();
	} else
	{
		std::cout << "Error writing output file for _array_statemonitor_1__indices." << endl;
	}
	ofstream outfile__array_statemonitor_1_clock_dt;
	outfile__array_statemonitor_1_clock_dt.open("results\\_array_statemonitor_1_clock_dt_3971040408683677175", ios::binary | ios::out);
	if(outfile__array_statemonitor_1_clock_dt.is_open())
	{
		outfile__array_statemonitor_1_clock_dt.write(reinterpret_cast<char*>(_array_statemonitor_1_clock_dt), 1*sizeof(_array_statemonitor_1_clock_dt[0]));
		outfile__array_statemonitor_1_clock_dt.close();
	} else
	{
		std::cout << "Error writing output file for _array_statemonitor_1_clock_dt." << endl;
	}
	ofstream outfile__array_statemonitor_1_clock_t;
	outfile__array_statemonitor_1_clock_t.open("results\\_array_statemonitor_1_clock_t_-7262499799254242934", ios::binary | ios::out);
	if(outfile__array_statemonitor_1_clock_t.is_open())
	{
		outfile__array_statemonitor_1_clock_t.write(reinterpret_cast<char*>(_array_statemonitor_1_clock_t), 1*sizeof(_array_statemonitor_1_clock_t[0]));
		outfile__array_statemonitor_1_clock_t.close();
	} else
	{
		std::cout << "Error writing output file for _array_statemonitor_1_clock_t." << endl;
	}
	ofstream outfile__array_statemonitor_1_clock_timestep;
	outfile__array_statemonitor_1_clock_timestep.open("results\\_array_statemonitor_1_clock_timestep_-6538729364457138717", ios::binary | ios::out);
	if(outfile__array_statemonitor_1_clock_timestep.is_open())
	{
		outfile__array_statemonitor_1_clock_timestep.write(reinterpret_cast<char*>(_array_statemonitor_1_clock_timestep), 1*sizeof(_array_statemonitor_1_clock_timestep[0]));
		outfile__array_statemonitor_1_clock_timestep.close();
	} else
	{
		std::cout << "Error writing output file for _array_statemonitor_1_clock_timestep." << endl;
	}
	ofstream outfile__array_statemonitor_1_N;
	outfile__array_statemonitor_1_N.open("results\\_array_statemonitor_1_N_-8106655974620711918", ios::binary | ios::out);
	if(outfile__array_statemonitor_1_N.is_open())
	{
		outfile__array_statemonitor_1_N.write(reinterpret_cast<char*>(_array_statemonitor_1_N), 1*sizeof(_array_statemonitor_1_N[0]));
		outfile__array_statemonitor_1_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_statemonitor_1_N." << endl;
	}
	ofstream outfile__array_statemonitor__indices;
	outfile__array_statemonitor__indices.open("results\\_array_statemonitor__indices_3017584569760417117", ios::binary | ios::out);
	if(outfile__array_statemonitor__indices.is_open())
	{
		outfile__array_statemonitor__indices.write(reinterpret_cast<char*>(_array_statemonitor__indices), 100*sizeof(_array_statemonitor__indices[0]));
		outfile__array_statemonitor__indices.close();
	} else
	{
		std::cout << "Error writing output file for _array_statemonitor__indices." << endl;
	}
	ofstream outfile__array_statemonitor_clock_dt;
	outfile__array_statemonitor_clock_dt.open("results\\_array_statemonitor_clock_dt_1386070905394900225", ios::binary | ios::out);
	if(outfile__array_statemonitor_clock_dt.is_open())
	{
		outfile__array_statemonitor_clock_dt.write(reinterpret_cast<char*>(_array_statemonitor_clock_dt), 1*sizeof(_array_statemonitor_clock_dt[0]));
		outfile__array_statemonitor_clock_dt.close();
	} else
	{
		std::cout << "Error writing output file for _array_statemonitor_clock_dt." << endl;
	}
	ofstream outfile__array_statemonitor_clock_t;
	outfile__array_statemonitor_clock_t.open("results\\_array_statemonitor_clock_t_5650692798298943133", ios::binary | ios::out);
	if(outfile__array_statemonitor_clock_t.is_open())
	{
		outfile__array_statemonitor_clock_t.write(reinterpret_cast<char*>(_array_statemonitor_clock_t), 1*sizeof(_array_statemonitor_clock_t[0]));
		outfile__array_statemonitor_clock_t.close();
	} else
	{
		std::cout << "Error writing output file for _array_statemonitor_clock_t." << endl;
	}
	ofstream outfile__array_statemonitor_clock_timestep;
	outfile__array_statemonitor_clock_timestep.open("results\\_array_statemonitor_clock_timestep_-6848382139332618038", ios::binary | ios::out);
	if(outfile__array_statemonitor_clock_timestep.is_open())
	{
		outfile__array_statemonitor_clock_timestep.write(reinterpret_cast<char*>(_array_statemonitor_clock_timestep), 1*sizeof(_array_statemonitor_clock_timestep[0]));
		outfile__array_statemonitor_clock_timestep.close();
	} else
	{
		std::cout << "Error writing output file for _array_statemonitor_clock_timestep." << endl;
	}
	ofstream outfile__array_statemonitor_N;
	outfile__array_statemonitor_N.open("results\\_array_statemonitor_N_-5645085762349242235", ios::binary | ios::out);
	if(outfile__array_statemonitor_N.is_open())
	{
		outfile__array_statemonitor_N.write(reinterpret_cast<char*>(_array_statemonitor_N), 1*sizeof(_array_statemonitor_N[0]));
		outfile__array_statemonitor_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_statemonitor_N." << endl;
	}
	ofstream outfile__array_synapses_N;
	outfile__array_synapses_N.open("results\\_array_synapses_N_4639461291556732226", ios::binary | ios::out);
	if(outfile__array_synapses_N.is_open())
	{
		outfile__array_synapses_N.write(reinterpret_cast<char*>(_array_synapses_N), 1*sizeof(_array_synapses_N[0]));
		outfile__array_synapses_N.close();
	} else
	{
		std::cout << "Error writing output file for _array_synapses_N." << endl;
	}

	ofstream outfile__dynamic_array_spikemonitor_1_i;
	outfile__dynamic_array_spikemonitor_1_i.open("results\\_dynamic_array_spikemonitor_1_i_2044251474265767121", ios::binary | ios::out);
	if(outfile__dynamic_array_spikemonitor_1_i.is_open())
	{
        if (! _dynamic_array_spikemonitor_1_i.empty() )
        {
			outfile__dynamic_array_spikemonitor_1_i.write(reinterpret_cast<char*>(&_dynamic_array_spikemonitor_1_i[0]), _dynamic_array_spikemonitor_1_i.size()*sizeof(_dynamic_array_spikemonitor_1_i[0]));
		    outfile__dynamic_array_spikemonitor_1_i.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_spikemonitor_1_i." << endl;
	}
	ofstream outfile__dynamic_array_spikemonitor_1_t;
	outfile__dynamic_array_spikemonitor_1_t.open("results\\_dynamic_array_spikemonitor_1_t_-13403564174063528", ios::binary | ios::out);
	if(outfile__dynamic_array_spikemonitor_1_t.is_open())
	{
        if (! _dynamic_array_spikemonitor_1_t.empty() )
        {
			outfile__dynamic_array_spikemonitor_1_t.write(reinterpret_cast<char*>(&_dynamic_array_spikemonitor_1_t[0]), _dynamic_array_spikemonitor_1_t.size()*sizeof(_dynamic_array_spikemonitor_1_t[0]));
		    outfile__dynamic_array_spikemonitor_1_t.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_spikemonitor_1_t." << endl;
	}
	ofstream outfile__dynamic_array_spikemonitor_i;
	outfile__dynamic_array_spikemonitor_i.open("results\\_dynamic_array_spikemonitor_i_6725483051249551373", ios::binary | ios::out);
	if(outfile__dynamic_array_spikemonitor_i.is_open())
	{
        if (! _dynamic_array_spikemonitor_i.empty() )
        {
			outfile__dynamic_array_spikemonitor_i.write(reinterpret_cast<char*>(&_dynamic_array_spikemonitor_i[0]), _dynamic_array_spikemonitor_i.size()*sizeof(_dynamic_array_spikemonitor_i[0]));
		    outfile__dynamic_array_spikemonitor_i.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_spikemonitor_i." << endl;
	}
	ofstream outfile__dynamic_array_spikemonitor_t;
	outfile__dynamic_array_spikemonitor_t.open("results\\_dynamic_array_spikemonitor_t_5166449136247538874", ios::binary | ios::out);
	if(outfile__dynamic_array_spikemonitor_t.is_open())
	{
        if (! _dynamic_array_spikemonitor_t.empty() )
        {
			outfile__dynamic_array_spikemonitor_t.write(reinterpret_cast<char*>(&_dynamic_array_spikemonitor_t[0]), _dynamic_array_spikemonitor_t.size()*sizeof(_dynamic_array_spikemonitor_t[0]));
		    outfile__dynamic_array_spikemonitor_t.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_spikemonitor_t." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_1_t;
	outfile__dynamic_array_statemonitor_1_t.open("results\\_dynamic_array_statemonitor_1_t_-6733914014273985203", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_1_t.is_open())
	{
        if (! _dynamic_array_statemonitor_1_t.empty() )
        {
			outfile__dynamic_array_statemonitor_1_t.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_1_t[0]), _dynamic_array_statemonitor_1_t.size()*sizeof(_dynamic_array_statemonitor_1_t[0]));
		    outfile__dynamic_array_statemonitor_1_t.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_1_t." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_t;
	outfile__dynamic_array_statemonitor_t.open("results\\_dynamic_array_statemonitor_t_-7456163220329060995", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_t.is_open())
	{
        if (! _dynamic_array_statemonitor_t.empty() )
        {
			outfile__dynamic_array_statemonitor_t.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_t[0]), _dynamic_array_statemonitor_t.size()*sizeof(_dynamic_array_statemonitor_t[0]));
		    outfile__dynamic_array_statemonitor_t.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_t." << endl;
	}
	ofstream outfile__dynamic_array_synapses__synaptic_post;
	outfile__dynamic_array_synapses__synaptic_post.open("results\\_dynamic_array_synapses__synaptic_post_5057333223235559440", ios::binary | ios::out);
	if(outfile__dynamic_array_synapses__synaptic_post.is_open())
	{
        if (! _dynamic_array_synapses__synaptic_post.empty() )
        {
			outfile__dynamic_array_synapses__synaptic_post.write(reinterpret_cast<char*>(&_dynamic_array_synapses__synaptic_post[0]), _dynamic_array_synapses__synaptic_post.size()*sizeof(_dynamic_array_synapses__synaptic_post[0]));
		    outfile__dynamic_array_synapses__synaptic_post.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_synapses__synaptic_post." << endl;
	}
	ofstream outfile__dynamic_array_synapses__synaptic_pre;
	outfile__dynamic_array_synapses__synaptic_pre.open("results\\_dynamic_array_synapses__synaptic_pre_5248430759710140046", ios::binary | ios::out);
	if(outfile__dynamic_array_synapses__synaptic_pre.is_open())
	{
        if (! _dynamic_array_synapses__synaptic_pre.empty() )
        {
			outfile__dynamic_array_synapses__synaptic_pre.write(reinterpret_cast<char*>(&_dynamic_array_synapses__synaptic_pre[0]), _dynamic_array_synapses__synaptic_pre.size()*sizeof(_dynamic_array_synapses__synaptic_pre[0]));
		    outfile__dynamic_array_synapses__synaptic_pre.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_synapses__synaptic_pre." << endl;
	}
	ofstream outfile__dynamic_array_synapses_apost;
	outfile__dynamic_array_synapses_apost.open("results\\_dynamic_array_synapses_apost_405275053876962530", ios::binary | ios::out);
	if(outfile__dynamic_array_synapses_apost.is_open())
	{
        if (! _dynamic_array_synapses_apost.empty() )
        {
			outfile__dynamic_array_synapses_apost.write(reinterpret_cast<char*>(&_dynamic_array_synapses_apost[0]), _dynamic_array_synapses_apost.size()*sizeof(_dynamic_array_synapses_apost[0]));
		    outfile__dynamic_array_synapses_apost.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_synapses_apost." << endl;
	}
	ofstream outfile__dynamic_array_synapses_apre;
	outfile__dynamic_array_synapses_apre.open("results\\_dynamic_array_synapses_apre_-7603486010187138516", ios::binary | ios::out);
	if(outfile__dynamic_array_synapses_apre.is_open())
	{
        if (! _dynamic_array_synapses_apre.empty() )
        {
			outfile__dynamic_array_synapses_apre.write(reinterpret_cast<char*>(&_dynamic_array_synapses_apre[0]), _dynamic_array_synapses_apre.size()*sizeof(_dynamic_array_synapses_apre[0]));
		    outfile__dynamic_array_synapses_apre.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_synapses_apre." << endl;
	}
	ofstream outfile__dynamic_array_synapses_delay;
	outfile__dynamic_array_synapses_delay.open("results\\_dynamic_array_synapses_delay_-3820318948673352934", ios::binary | ios::out);
	if(outfile__dynamic_array_synapses_delay.is_open())
	{
        if (! _dynamic_array_synapses_delay.empty() )
        {
			outfile__dynamic_array_synapses_delay.write(reinterpret_cast<char*>(&_dynamic_array_synapses_delay[0]), _dynamic_array_synapses_delay.size()*sizeof(_dynamic_array_synapses_delay[0]));
		    outfile__dynamic_array_synapses_delay.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_synapses_delay." << endl;
	}
	ofstream outfile__dynamic_array_synapses_delay_1;
	outfile__dynamic_array_synapses_delay_1.open("results\\_dynamic_array_synapses_delay_1_6636737312748171074", ios::binary | ios::out);
	if(outfile__dynamic_array_synapses_delay_1.is_open())
	{
        if (! _dynamic_array_synapses_delay_1.empty() )
        {
			outfile__dynamic_array_synapses_delay_1.write(reinterpret_cast<char*>(&_dynamic_array_synapses_delay_1[0]), _dynamic_array_synapses_delay_1.size()*sizeof(_dynamic_array_synapses_delay_1[0]));
		    outfile__dynamic_array_synapses_delay_1.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_synapses_delay_1." << endl;
	}
	ofstream outfile__dynamic_array_synapses_lastupdate;
	outfile__dynamic_array_synapses_lastupdate.open("results\\_dynamic_array_synapses_lastupdate_-2721785242155238857", ios::binary | ios::out);
	if(outfile__dynamic_array_synapses_lastupdate.is_open())
	{
        if (! _dynamic_array_synapses_lastupdate.empty() )
        {
			outfile__dynamic_array_synapses_lastupdate.write(reinterpret_cast<char*>(&_dynamic_array_synapses_lastupdate[0]), _dynamic_array_synapses_lastupdate.size()*sizeof(_dynamic_array_synapses_lastupdate[0]));
		    outfile__dynamic_array_synapses_lastupdate.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_synapses_lastupdate." << endl;
	}
	ofstream outfile__dynamic_array_synapses_N_incoming;
	outfile__dynamic_array_synapses_N_incoming.open("results\\_dynamic_array_synapses_N_incoming_-726932618916879662", ios::binary | ios::out);
	if(outfile__dynamic_array_synapses_N_incoming.is_open())
	{
        if (! _dynamic_array_synapses_N_incoming.empty() )
        {
			outfile__dynamic_array_synapses_N_incoming.write(reinterpret_cast<char*>(&_dynamic_array_synapses_N_incoming[0]), _dynamic_array_synapses_N_incoming.size()*sizeof(_dynamic_array_synapses_N_incoming[0]));
		    outfile__dynamic_array_synapses_N_incoming.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_synapses_N_incoming." << endl;
	}
	ofstream outfile__dynamic_array_synapses_N_outgoing;
	outfile__dynamic_array_synapses_N_outgoing.open("results\\_dynamic_array_synapses_N_outgoing_-5272174316927614017", ios::binary | ios::out);
	if(outfile__dynamic_array_synapses_N_outgoing.is_open())
	{
        if (! _dynamic_array_synapses_N_outgoing.empty() )
        {
			outfile__dynamic_array_synapses_N_outgoing.write(reinterpret_cast<char*>(&_dynamic_array_synapses_N_outgoing[0]), _dynamic_array_synapses_N_outgoing.size()*sizeof(_dynamic_array_synapses_N_outgoing[0]));
		    outfile__dynamic_array_synapses_N_outgoing.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_synapses_N_outgoing." << endl;
	}
	ofstream outfile__dynamic_array_synapses_R_hat;
	outfile__dynamic_array_synapses_R_hat.open("results\\_dynamic_array_synapses_R_hat_4399214112561688751", ios::binary | ios::out);
	if(outfile__dynamic_array_synapses_R_hat.is_open())
	{
        if (! _dynamic_array_synapses_R_hat.empty() )
        {
			outfile__dynamic_array_synapses_R_hat.write(reinterpret_cast<char*>(&_dynamic_array_synapses_R_hat[0]), _dynamic_array_synapses_R_hat.size()*sizeof(_dynamic_array_synapses_R_hat[0]));
		    outfile__dynamic_array_synapses_R_hat.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_synapses_R_hat." << endl;
	}
	ofstream outfile__dynamic_array_synapses_w;
	outfile__dynamic_array_synapses_w.open("results\\_dynamic_array_synapses_w_8513689814937615011", ios::binary | ios::out);
	if(outfile__dynamic_array_synapses_w.is_open())
	{
        if (! _dynamic_array_synapses_w.empty() )
        {
			outfile__dynamic_array_synapses_w.write(reinterpret_cast<char*>(&_dynamic_array_synapses_w[0]), _dynamic_array_synapses_w.size()*sizeof(_dynamic_array_synapses_w[0]));
		    outfile__dynamic_array_synapses_w.close();
		}
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_synapses_w." << endl;
	}

	ofstream outfile__dynamic_array_statemonitor_1_L;
	outfile__dynamic_array_statemonitor_1_L.open("results\\_dynamic_array_statemonitor_1_L_6079773891958699937", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_1_L.is_open())
	{
        for (int n=0; n<_dynamic_array_statemonitor_1_L.n; n++)
        {
            if (! _dynamic_array_statemonitor_1_L(n).empty())
            {
                outfile__dynamic_array_statemonitor_1_L.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_1_L(n, 0)), _dynamic_array_statemonitor_1_L.m*sizeof(_dynamic_array_statemonitor_1_L(0, 0)));
            }
        }
        outfile__dynamic_array_statemonitor_1_L.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_1_L." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_1_u;
	outfile__dynamic_array_statemonitor_1_u.open("results\\_dynamic_array_statemonitor_1_u_-457272687482333811", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_1_u.is_open())
	{
        for (int n=0; n<_dynamic_array_statemonitor_1_u.n; n++)
        {
            if (! _dynamic_array_statemonitor_1_u(n).empty())
            {
                outfile__dynamic_array_statemonitor_1_u.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_1_u(n, 0)), _dynamic_array_statemonitor_1_u.m*sizeof(_dynamic_array_statemonitor_1_u(0, 0)));
            }
        }
        outfile__dynamic_array_statemonitor_1_u.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_1_u." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_1_v;
	outfile__dynamic_array_statemonitor_1_v.open("results\\_dynamic_array_statemonitor_1_v_4439588890072526155", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_1_v.is_open())
	{
        for (int n=0; n<_dynamic_array_statemonitor_1_v.n; n++)
        {
            if (! _dynamic_array_statemonitor_1_v(n).empty())
            {
                outfile__dynamic_array_statemonitor_1_v.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_1_v(n, 0)), _dynamic_array_statemonitor_1_v.m*sizeof(_dynamic_array_statemonitor_1_v(0, 0)));
            }
        }
        outfile__dynamic_array_statemonitor_1_v.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_1_v." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_K;
	outfile__dynamic_array_statemonitor_K.open("results\\_dynamic_array_statemonitor_K_7861085328579309433", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_K.is_open())
	{
        for (int n=0; n<_dynamic_array_statemonitor_K.n; n++)
        {
            if (! _dynamic_array_statemonitor_K(n).empty())
            {
                outfile__dynamic_array_statemonitor_K.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_K(n, 0)), _dynamic_array_statemonitor_K.m*sizeof(_dynamic_array_statemonitor_K(0, 0)));
            }
        }
        outfile__dynamic_array_statemonitor_K.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_K." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_R_hat;
	outfile__dynamic_array_statemonitor_R_hat.open("results\\_dynamic_array_statemonitor_R_hat_3290712765598294941", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_R_hat.is_open())
	{
        for (int n=0; n<_dynamic_array_statemonitor_R_hat.n; n++)
        {
            if (! _dynamic_array_statemonitor_R_hat(n).empty())
            {
                outfile__dynamic_array_statemonitor_R_hat.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_R_hat(n, 0)), _dynamic_array_statemonitor_R_hat.m*sizeof(_dynamic_array_statemonitor_R_hat(0, 0)));
            }
        }
        outfile__dynamic_array_statemonitor_R_hat.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_R_hat." << endl;
	}
	ofstream outfile__dynamic_array_statemonitor_w;
	outfile__dynamic_array_statemonitor_w.open("results\\_dynamic_array_statemonitor_w_2540624823645533936", ios::binary | ios::out);
	if(outfile__dynamic_array_statemonitor_w.is_open())
	{
        for (int n=0; n<_dynamic_array_statemonitor_w.n; n++)
        {
            if (! _dynamic_array_statemonitor_w(n).empty())
            {
                outfile__dynamic_array_statemonitor_w.write(reinterpret_cast<char*>(&_dynamic_array_statemonitor_w(n, 0)), _dynamic_array_statemonitor_w.m*sizeof(_dynamic_array_statemonitor_w(0, 0)));
            }
        }
        outfile__dynamic_array_statemonitor_w.close();
	} else
	{
		std::cout << "Error writing output file for _dynamic_array_statemonitor_w." << endl;
	}

	// Write profiling info to disk
	ofstream outfile_profiling_info;
	outfile_profiling_info.open("results/profiling_info.txt", ios::out);
	if(outfile_profiling_info.is_open())
	{
	outfile_profiling_info << "neurongroup_resetter_codeobject\t" << neurongroup_resetter_codeobject_profiling_info << std::endl;
	outfile_profiling_info << "neurongroup_stateupdater_codeobject\t" << neurongroup_stateupdater_codeobject_profiling_info << std::endl;
	outfile_profiling_info << "neurongroup_thresholder_codeobject\t" << neurongroup_thresholder_codeobject_profiling_info << std::endl;
	outfile_profiling_info << "poissongroup_thresholder_codeobject\t" << poissongroup_thresholder_codeobject_profiling_info << std::endl;
	outfile_profiling_info << "spikemonitor_1_codeobject\t" << spikemonitor_1_codeobject_profiling_info << std::endl;
	outfile_profiling_info << "spikemonitor_codeobject\t" << spikemonitor_codeobject_profiling_info << std::endl;
	outfile_profiling_info << "statemonitor_1_codeobject\t" << statemonitor_1_codeobject_profiling_info << std::endl;
	outfile_profiling_info << "statemonitor_codeobject\t" << statemonitor_codeobject_profiling_info << std::endl;
	outfile_profiling_info << "synapses_post_codeobject\t" << synapses_post_codeobject_profiling_info << std::endl;
	outfile_profiling_info << "synapses_post_initialise_queue\t" << synapses_post_initialise_queue_profiling_info << std::endl;
	outfile_profiling_info << "synapses_post_push_spikes\t" << synapses_post_push_spikes_profiling_info << std::endl;
	outfile_profiling_info << "synapses_pre_codeobject\t" << synapses_pre_codeobject_profiling_info << std::endl;
	outfile_profiling_info << "synapses_pre_initialise_queue\t" << synapses_pre_initialise_queue_profiling_info << std::endl;
	outfile_profiling_info << "synapses_pre_push_spikes\t" << synapses_pre_push_spikes_profiling_info << std::endl;
	outfile_profiling_info << "synapses_stateupdater_codeobject\t" << synapses_stateupdater_codeobject_profiling_info << std::endl;
	outfile_profiling_info << "synapses_synapses_create_generator_codeobject\t" << synapses_synapses_create_generator_codeobject_profiling_info << std::endl;
	outfile_profiling_info.close();
	} else
	{
	    std::cout << "Error writing profiling info to file." << std::endl;
	}

	// Write last run info to disk
	ofstream outfile_last_run_info;
	outfile_last_run_info.open("results/last_run_info.txt", ios::out);
	if(outfile_last_run_info.is_open())
	{
		outfile_last_run_info << (Network::_last_run_time) << " " << (Network::_last_run_completed_fraction) << std::endl;
		outfile_last_run_info.close();
	} else
	{
	    std::cout << "Error writing last run info to file." << std::endl;
	}
}

void _dealloc_arrays()
{
	using namespace brian;


	// static arrays
	if(_static_array__array_poissongroup_rates!=0)
	{
		delete [] _static_array__array_poissongroup_rates;
		_static_array__array_poissongroup_rates = 0;
	}
	if(_static_array__array_statemonitor__indices!=0)
	{
		delete [] _static_array__array_statemonitor__indices;
		_static_array__array_statemonitor__indices = 0;
	}
}

