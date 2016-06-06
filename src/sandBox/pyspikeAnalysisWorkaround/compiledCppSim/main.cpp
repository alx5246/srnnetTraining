#include <stdlib.h>
#include "objects.h"
#include <ctime>
#include <time.h>

#include "run.h"
#include "brianlib/common_math.h"

#include "code_objects/synapses_post_initialise_queue.h"
#include "code_objects/synapses_stateupdater_codeobject.h"
#include "code_objects/neurongroup_thresholder_codeobject.h"
#include "code_objects/spikemonitor_codeobject.h"
#include "code_objects/synapses_post_codeobject.h"
#include "code_objects/neurongroup_resetter_codeobject.h"
#include "code_objects/synapses_pre_initialise_queue.h"
#include "code_objects/statemonitor_1_codeobject.h"
#include "code_objects/synapses_group_variable_set_conditional_codeobject.h"
#include "code_objects/poissongroup_thresholder_codeobject.h"
#include "code_objects/synapses_synapses_create_generator_codeobject.h"
#include "code_objects/synapses_pre_push_spikes.h"
#include "code_objects/neurongroup_stateupdater_codeobject.h"
#include "code_objects/statemonitor_codeobject.h"
#include "code_objects/spikemonitor_1_codeobject.h"
#include "code_objects/synapses_pre_codeobject.h"
#include "code_objects/synapses_post_push_spikes.h"


#include <iostream>
#include <fstream>




int main(int argc, char **argv)
{

	brian_start();

	{
		using namespace brian;

		
                
        _array_defaultclock_dt[0] = 0.0001;
        _array_neurongroup_lastspike[0] = - INFINITY;
        _array_neurongroup_not_refractory[0] = true;
        
                        
                        for(int i=0; i<_num__array_poissongroup_rates; i++)
                        {
                            _array_poissongroup_rates[i] = _static_array__array_poissongroup_rates[i];
                        }
                        
        _run_synapses_synapses_create_generator_codeobject();
        
                        
                        for(int i=0; i<_dynamic_array_synapses_w.size(); i++)
                        {
                            _dynamic_array_synapses_w[i] = 2.0;
                        }
                        
        _array_statemonitor_clock_dt[0] = 0.0005;
        
                        
                        for(int i=0; i<_num__array_statemonitor__indices; i++)
                        {
                            _array_statemonitor__indices[i] = _static_array__array_statemonitor__indices[i];
                        }
                        
        _array_statemonitor_1_clock_dt[0] = 0.0005;
        _array_statemonitor_1__indices[0] = 0;
        _array_statemonitor_1_clock_timestep[0] = 0;
        _array_statemonitor_1_clock_t[0] = 0.0;
        _array_statemonitor_clock_timestep[0] = 0;
        _array_statemonitor_clock_t[0] = 0.0;
        _array_defaultclock_dt[0] = 9.999999999999999e-05;
        _array_defaultclock_timestep[0] = 0;
        _array_defaultclock_t[0] = 0.0;
        _run_synapses_group_variable_set_conditional_codeobject();
        _run_synapses_pre_initialise_queue();
        _run_synapses_post_initialise_queue();
        magicnetwork.clear();
        magicnetwork.add(&statemonitor_clock, _run_statemonitor_codeobject);
        magicnetwork.add(&statemonitor_1_clock, _run_statemonitor_1_codeobject);
        magicnetwork.add(&defaultclock, _run_neurongroup_stateupdater_codeobject);
        magicnetwork.add(&defaultclock, _run_synapses_stateupdater_codeobject);
        magicnetwork.add(&defaultclock, _run_neurongroup_thresholder_codeobject);
        magicnetwork.add(&defaultclock, _run_poissongroup_thresholder_codeobject);
        magicnetwork.add(&defaultclock, _run_spikemonitor_codeobject);
        magicnetwork.add(&defaultclock, _run_spikemonitor_1_codeobject);
        magicnetwork.add(&defaultclock, _run_synapses_pre_push_spikes);
        magicnetwork.add(&defaultclock, _run_synapses_pre_codeobject);
        magicnetwork.add(&defaultclock, _run_synapses_post_push_spikes);
        magicnetwork.add(&defaultclock, _run_synapses_post_codeobject);
        magicnetwork.add(&defaultclock, _run_neurongroup_resetter_codeobject);
        magicnetwork.run(1.0, NULL, 10.0);
        _debugmsg_spikemonitor_codeobject();
        
        _debugmsg_synapses_post_codeobject();
        
        _debugmsg_spikemonitor_1_codeobject();
        
        _debugmsg_synapses_pre_codeobject();

	}

	brian_end();

	return 0;
}