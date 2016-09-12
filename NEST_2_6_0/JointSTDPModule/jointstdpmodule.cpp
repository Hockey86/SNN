/*
 *  jointstdpmodule.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

// include necessary NEST headers
#include "config.h"
#include "network.h"
#include "model.h"
#include "dynamicloader.h"
#include "genericmodel.h"
#include "booldatum.h"
#include "integerdatum.h"
#include "tokenarray.h"
#include "exceptions.h"
#include "sliexceptions.h"
#include "nestmodule.h"
#include "connector_model_impl.h"
#include "target_identifier.h"

// include headers with your own stuff
#include "jointstdpmodule.h"
#include "parrot_neuron2.h"
#include "spike_response_neuron.h"
//#include "weight_stdp_axonal_connection.h"
//#include "joint_stdp_axonal_connection.h"
#include "weight_stdp_connection.h"
#include "joint_stdp_connection.h"
#include "immediate_inhibition_connection.h"

// -- Interface to dynamic module loader ---------------------------------------

/*
 * The dynamic module loader must be able to find your module.
 * You make the module known to the loader by defining an instance of your
 * module class in global scope. This instance must have the name
 *
 * <modulename>_LTX_mod
 *
 * The dynamicloader can then load modulename and search for symbol "mod" in it.
 */

mynest::JointSTDPModule jointstdpmodule_LTX_mod;

// -- DynModule functions ------------------------------------------------------

mynest::JointSTDPModule::JointSTDPModule()
  {
#ifdef LINKED_MODULE
     // register this module at the dynamic loader
     // this is needed to allow for linking in this module at compile time
     // all registered modules will be initialized by the main app's dynamic loader
     nest::DynamicLoaderModule::registerLinkedModule(this);
#endif
   }

mynest::JointSTDPModule::~JointSTDPModule()
   {}

   const std::string mynest::JointSTDPModule::name(void) const
   {
     return std::string("Joint STDP Module"); // Return name of the module
   }

   const std::string mynest::JointSTDPModule::commandstring(void) const
   {
     // Instruct the interpreter to load jointstdpmodule-init.sli
     return std::string("(jointstdpmodule-init) run");
   }

  //-------------------------------------------------------------------------------------

  void mynest::JointSTDPModule::init(SLIInterpreter *i, nest::Network*)
  {
    /* Register a neuron or device model.
       Give node type as template argument and the name as second argument.
       The first argument is always a reference to the network.
       Return value is a handle for later unregistration.
    */
    nest::register_model<parrot_neuron2>(nest::NestModule::get_network(),
                                        "parrot_neuron2");
    nest::register_model<spike_response_neuron>(nest::NestModule::get_network(),
                                        "spike_response_neuron");

    /* Register a synapse type.
       Give synapse type as template argument and the name as second argument.
       The first argument is always a reference to the network.

       There are two choices for the template argument:
       	   - nest::TargetIdentifierPtrRport
       	   - nest::TargetIdentifierIndex
       The first is the standard and you should usually stick to it.
       nest::TargetIdentifierIndex reduces the memory requirement of synapses
       even further, but limits the number of available rports. Please see
       Kunkel et al, Front Neurofinfom 8:78 (2014), Sec 3.3.2, for details.
    */
    //nest::register_connection_model<WeightSTDPAxonalConnection<nest::TargetIdentifierPtrRport> >(nest::NestModule::get_network(),
		//								       "weight_stdp_axonal_synapse");
    //nest::register_connection_model<JointSTDPAxonalConnection<nest::TargetIdentifierPtrRport> >(nest::NestModule::get_network(),
		//								       "joint_stdp_axonal_synapse");
    nest::register_connection_model<WeightSTDPConnection<nest::TargetIdentifierPtrRport> >(nest::NestModule::get_network(),
										       "weight_stdp_synapse");
    nest::register_connection_model<JointSTDPConnection<nest::TargetIdentifierPtrRport> >(nest::NestModule::get_network(),
										       "joint_stdp_synapse");
    nest::register_connection_model<ImmediateInhibitionConnection<nest::TargetIdentifierPtrRport> >(nest::NestModule::get_network(),
										       "immediate_inhibition_synapse");

    /* Register a SLI function.
       The first argument is the function name for SLI, the second a pointer to
       the function object. If you do not want to overload the function in SLI,
       you do not need to give the mangled name. If you give a mangled name, you
       should define a type trie in the jointstdpmodule-init.sli file.
    */
    //i->createcommand("StepPatternConnect_Vi_i_Vi_i_l",&stepPatternConnect_Vi_i_Vi_i_lFunction);

  }  // JointSTDPModule::init()


