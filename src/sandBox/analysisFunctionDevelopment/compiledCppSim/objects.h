
#ifndef _BRIAN_OBJECTS_H
#define _BRIAN_OBJECTS_H

#include "synapses_classes.h"
#include "brianlib/clocks.h"
#include "brianlib/dynamic_array.h"
#include "brianlib/stdint_compat.h"
#include "network.h"
#include "randomkit.h"
#include<vector>


namespace brian {

// In OpenMP we need one state per thread
extern std::vector< rk_state* > _mersenne_twister_states;

//////////////// clocks ///////////////////
extern Clock defaultclock;
extern Clock statemonitor_clock;
extern Clock statemonitor_1_clock;

//////////////// networks /////////////////
extern Network magicnetwork;
extern Network magicnetwork;

//////////////// dynamic arrays ///////////
extern std::vector<int32_t> _dynamic_array_spikemonitor_1_i;
extern std::vector<double> _dynamic_array_spikemonitor_1_t;
extern std::vector<int32_t> _dynamic_array_spikemonitor_i;
extern std::vector<double> _dynamic_array_spikemonitor_t;
extern std::vector<double> _dynamic_array_statemonitor_1_t;
extern std::vector<double> _dynamic_array_statemonitor_t;
extern std::vector<int32_t> _dynamic_array_synapses__synaptic_post;
extern std::vector<int32_t> _dynamic_array_synapses__synaptic_pre;
extern std::vector<double> _dynamic_array_synapses_apost;
extern std::vector<double> _dynamic_array_synapses_apre;
extern std::vector<double> _dynamic_array_synapses_delay;
extern std::vector<double> _dynamic_array_synapses_delay_1;
extern std::vector<double> _dynamic_array_synapses_lastupdate;
extern std::vector<int32_t> _dynamic_array_synapses_N_incoming;
extern std::vector<int32_t> _dynamic_array_synapses_N_outgoing;
extern std::vector<double> _dynamic_array_synapses_R_hat;
extern std::vector<double> _dynamic_array_synapses_w;

//////////////// arrays ///////////////////
extern double *_array_defaultclock_dt;
extern const int _num__array_defaultclock_dt;
extern double *_array_defaultclock_t;
extern const int _num__array_defaultclock_t;
extern uint64_t *_array_defaultclock_timestep;
extern const int _num__array_defaultclock_timestep;
extern int32_t *_array_neurongroup__spikespace;
extern const int _num__array_neurongroup__spikespace;
extern int32_t *_array_neurongroup_i;
extern const int _num__array_neurongroup_i;
extern double *_array_neurongroup_I;
extern const int _num__array_neurongroup_I;
extern double *_array_neurongroup_L;
extern const int _num__array_neurongroup_L;
extern double *_array_neurongroup_lastspike;
extern const int _num__array_neurongroup_lastspike;
extern bool *_array_neurongroup_not_refractory;
extern const int _num__array_neurongroup_not_refractory;
extern double *_array_neurongroup_u;
extern const int _num__array_neurongroup_u;
extern double *_array_neurongroup_v;
extern const int _num__array_neurongroup_v;
extern int32_t *_array_poissongroup__spikespace;
extern const int _num__array_poissongroup__spikespace;
extern int32_t *_array_poissongroup_i;
extern const int _num__array_poissongroup_i;
extern double *_array_poissongroup_rates;
extern const int _num__array_poissongroup_rates;
extern int32_t *_array_spikemonitor_1__source_idx;
extern const int _num__array_spikemonitor_1__source_idx;
extern int32_t *_array_spikemonitor_1_count;
extern const int _num__array_spikemonitor_1_count;
extern int32_t *_array_spikemonitor_1_N;
extern const int _num__array_spikemonitor_1_N;
extern int32_t *_array_spikemonitor__source_idx;
extern const int _num__array_spikemonitor__source_idx;
extern int32_t *_array_spikemonitor_count;
extern const int _num__array_spikemonitor_count;
extern int32_t *_array_spikemonitor_N;
extern const int _num__array_spikemonitor_N;
extern int32_t *_array_statemonitor_1__indices;
extern const int _num__array_statemonitor_1__indices;
extern double *_array_statemonitor_1_clock_dt;
extern const int _num__array_statemonitor_1_clock_dt;
extern double *_array_statemonitor_1_clock_t;
extern const int _num__array_statemonitor_1_clock_t;
extern uint64_t *_array_statemonitor_1_clock_timestep;
extern const int _num__array_statemonitor_1_clock_timestep;
extern double *_array_statemonitor_1_L;
extern const int _num__array_statemonitor_1_L;
extern int32_t *_array_statemonitor_1_N;
extern const int _num__array_statemonitor_1_N;
extern double *_array_statemonitor_1_u;
extern const int _num__array_statemonitor_1_u;
extern double *_array_statemonitor_1_v;
extern const int _num__array_statemonitor_1_v;
extern int32_t *_array_statemonitor__indices;
extern const int _num__array_statemonitor__indices;
extern double *_array_statemonitor_clock_dt;
extern const int _num__array_statemonitor_clock_dt;
extern double *_array_statemonitor_clock_t;
extern const int _num__array_statemonitor_clock_t;
extern uint64_t *_array_statemonitor_clock_timestep;
extern const int _num__array_statemonitor_clock_timestep;
extern double *_array_statemonitor_K;
extern const int _num__array_statemonitor_K;
extern int32_t *_array_statemonitor_N;
extern const int _num__array_statemonitor_N;
extern double *_array_statemonitor_R_hat;
extern const int _num__array_statemonitor_R_hat;
extern double *_array_statemonitor_w;
extern const int _num__array_statemonitor_w;
extern int32_t *_array_synapses_N;
extern const int _num__array_synapses_N;

//////////////// dynamic arrays 2d /////////
extern DynamicArray2D<double> _dynamic_array_statemonitor_1_L;
extern DynamicArray2D<double> _dynamic_array_statemonitor_1_u;
extern DynamicArray2D<double> _dynamic_array_statemonitor_1_v;
extern DynamicArray2D<double> _dynamic_array_statemonitor_K;
extern DynamicArray2D<double> _dynamic_array_statemonitor_R_hat;
extern DynamicArray2D<double> _dynamic_array_statemonitor_w;

/////////////// static arrays /////////////
extern double *_static_array__array_poissongroup_rates;
extern const int _num__static_array__array_poissongroup_rates;
extern int32_t *_static_array__array_statemonitor__indices;
extern const int _num__static_array__array_statemonitor__indices;

//////////////// synapses /////////////////
// synapses
extern SynapticPathway<double> synapses_post;
extern SynapticPathway<double> synapses_pre;

// Profiling information for each code object
extern double neurongroup_resetter_codeobject_profiling_info;
extern double neurongroup_stateupdater_codeobject_profiling_info;
extern double neurongroup_thresholder_codeobject_profiling_info;
extern double poissongroup_thresholder_codeobject_profiling_info;
extern double spikemonitor_1_codeobject_profiling_info;
extern double spikemonitor_codeobject_profiling_info;
extern double statemonitor_1_codeobject_profiling_info;
extern double statemonitor_codeobject_profiling_info;
extern double synapses_post_codeobject_profiling_info;
extern double synapses_post_initialise_queue_profiling_info;
extern double synapses_post_push_spikes_profiling_info;
extern double synapses_pre_codeobject_profiling_info;
extern double synapses_pre_initialise_queue_profiling_info;
extern double synapses_pre_push_spikes_profiling_info;
extern double synapses_stateupdater_codeobject_profiling_info;
extern double synapses_synapses_create_generator_codeobject_profiling_info;

}

void _init_arrays();
void _load_arrays();
void _write_arrays();
void _dealloc_arrays();

#endif


