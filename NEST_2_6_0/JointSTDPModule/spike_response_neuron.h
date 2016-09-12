/*
 *  spike_response_neuron.h
 *
 *  This file is part of JointSTDPModule.
 *
 */

#ifndef SPIKE_RESPONSE_NEURON_H
#define SPIKE_RESPONSE_NEURON_H

#include "nest.h"
#include "event.h"
#include "my_base_neuron.h"
#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"
#include "recordables_map.h"
#include "normal_randomdev.h"

namespace mynest
{
  //class nest::Network;

  /* BeginDocumentation
     Name: spike_response_neuron - Spike-response model used in Carandini et al 2007,
     inherited from my_base_neuron.

     Description:
     The membrane potential is the sum of stereotyped events: the postsynaptic
     potentials (V_syn), waveforms that include a spike and the subsequent
     after-hyperpolarization (V_spike) and Gaussian-distributed white noise.
  
     Inherited from my_base_neuron, it stores its upstream connections
     and downstream neurons with immediate inhibition connection.

     The postsynaptic potential is described by alpha function where where
     U_epsp is the maximal amplitude of the EPSP, tau_epsp is the time to
     peak of the EPSP and tau_ipsp is the time to peak of the IPSP.

     The spike waveform is described as a delta peak followed by a membrane
     potential reset and exponential decay. U_reset is the magnitude of the
     reset/after-hyperpolarization and tau_reset is the time constant of
     recovery from this hyperpolarization.

     The linear subthresold dynamics is integrated by the Exact
     Integration scheme [1]. The neuron dynamics is solved on the time
     grid given by the computation step size. Incoming as well as emitted
     spikes are forced to that grid.

     Note:
     The way the noise term was implemented in the original model makes it
     unsuitable for simulation in NEST. The workaround was to prepare the
     noise signal externally prior to simulation. The noise signal,
     if present, has to be at least as long as the simulation.

     Parameters: 
     The following parameters can be set in the status dictionary.

     t_ref          double - Duration of refractory period in ms. 
     tau_epsp       double - Membrane time constant in ms.
     tau_ipsp       double - Membrane time constant in ms.
     tau_reset      double - Refractory time constant in ms.
     U_epsp         double - Maximum amplitude of the EPSP. Normalized.
     U_ipsp         double - Maximum amplitude of the IPSP. Normalized.
     U_reset        double - Reset value of the membrane potential. Normalized.
     U_noise        double - Noise scale. Normalized.
     noise   vector<double>- Noise signal.
 
     References:
     [1] Carandini M, Horton JC, Sincich LC (2007) Thalamic filtering of retinal
     spike trains by postsynaptic summation. J Vis 7(14):20,1-11.
     [2] Rotter S & Diesmann M (1999) Exact simulation of time-invariant linear
     systems with applications to neuronal modeling. Biologial Cybernetics
     81:381-402.

     Sends: SpikeEvent

     Receives: SpikeEvent, DataLoggingRequest

     FirstVersion: May 2012
     Author: Thomas Heiberg, Birgit Kriener
     
     Modified: Nov 2015, Sun, Haoqi
  */

  /**
   * Neuron model used in Carandini et al 2007.
   */
  class spike_response_neuron: public my_base_neuron
  {
    
  public:        
    
    spike_response_neuron();
    spike_response_neuron(const spike_response_neuron&);

    /**
     * Import sets of overloaded virtual functions.
     * @see Technical Issues / Virtual Functions: Overriding, Overloading, and Hiding
     */
    using nest::Node::handle;
    using nest::Node::handles_test_event;

    nest::port send_test_event(nest::Node&, nest::rport, nest::synindex, bool);
    
    void handle(nest::SpikeEvent &);
    void handle(nest::DataLoggingRequest &);
    
    nest::port handles_test_event(nest::SpikeEvent &, nest::rport);
    nest::port handles_test_event(nest::DataLoggingRequest &, nest::rport);

    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);
    
    void immediate_inhibit(nest::double_t w)
    {
      if (S_.r_==0 && V_.ImmediateInhibitionCounts_>0) S_.r_ = V_.ImmediateInhibitionCounts_;
      S_.V_spike_ += w*P_.U_ipsp_;
    }
    
    void immediate_excite(nest::double_t w)
    { S_.V_input_ += w*P_.U_epsp_; }

  private:

    void init_node_(const nest::Node& proto);
    void init_state_(const nest::Node& proto);
    void init_buffers_();
    void calibrate();

    void update(const nest::Time &, const nest::long_t, const nest::long_t);

    // The next two classes need to be friends to access the State_ class/member
    friend class nest::RecordablesMap<spike_response_neuron>;
    friend class nest::UniversalDataLogger<spike_response_neuron>;

    // ----------------------------------------------------------------

    /**
     * State variables of the model.
     */
    struct State_
    {
      // state variables
      nest::double_t i_syn_ex_;  // postsynaptic current for exc. inputs, variable 1
      nest::double_t i_syn_in_;  // postsynaptic current for inh. inputs, variable 1
      nest::double_t V_syn_;	   // psp waveform, variable 2
      nest::double_t V_spike_;	 // post spike reset waveform, variable 3
      nest::double_t V_m_;       // membrane potential, variable 4
      nest::double_t V_input_;   // the pulse from immediate_excitation_synapse
      nest::int_t    r_;  //!< Number of refractory steps remainingS.

      nest::ulong_t position_;

      State_();  //!< Default initialization

      void get(DictionaryDatum &) const;
      void set(DictionaryDatum const &);
    };

    // ----------------------------------------------------------------

    /** 
     * Independent parameters of the model. 
     */
    struct Parameters_
    {    
      /** Refractory period in ms. */
      nest::double_t TauR_;
  
      /** Membrane time constant in ms. */
      nest::double_t tau_epsp_;
  
      /** Membrane time constant in ms. */
      nest::double_t tau_ipsp_;
  
      /** Period of being completely inhibited due to immediate_inhibition_synapse, in ms. */
      nest::double_t immediate_inhibition_duration_;
  
      /** whether to flush PSP after firing. */
      bool flush_PSP_after_firing_;

      /** Refractory time constant in ms. */
      nest::double_t tau_reset_;

      /** Resting potential. Normalized = 0.0. */
      nest::double_t E_L_;

      /** Threshold. Normalized = 1.0. */
      nest::double_t U_th_;

      /** Normalized maximum amplitude of the EPSP. */
      nest::double_t U_epsp_;

      /** Normalized maximum amplitude of the IPSP. */
      nest::double_t U_ipsp_;

      /** Normalized magnitude of the membrane potential reset. */
      nest::double_t U_reset_;

      /** Membrane capacitance. Note: Does not have any function currently. */
      nest::double_t C_;

      /** Noise scale. */
      nest::double_t U_noise_;

      /** Noise signal. */
      std::vector<double_t> noise_;

      Parameters_();  //!< Sets default parameter values

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary

      /** Set values from dictionary.
       * @returns Change in reversal potential E_L, to be passed to State_::set()
       * @note State is passed so that the position can be reset if the
       *       noise_ vector has been filled with new data.
       */
      void set(const DictionaryDatum&, State_& s);
    };
    

    // ----------------------------------------------------------------

    /**
     * Buffers of the model.
     */
    struct Buffers_
    {
      Buffers_(spike_response_neuron &);
      Buffers_(const Buffers_ &, spike_response_neuron &);

      /** buffers and sums up incoming spikes/currents */
      nest::RingBuffer spikes_ex_;
      nest::RingBuffer spikes_in_;
      nest::RingBuffer currents_;  

      //! Logger for all analog data
      nest::UniversalDataLogger<spike_response_neuron> logger_;
    };
    
    // ---------------------------------------------------------------- 

    /**
     * Internal variables of the model.
     */
    struct Variables_
    { 
      /** Amplitude of the synaptic current.
	  This value is chosen such that a post-synaptic potential with
	  weight one has an amplitude of 1 mV.
	  @note mog - I assume this, not checked. 
      */
      //    double_t PSCInitialValue_;
    
      // time evolution operator
      nest::double_t P20_;
      nest::double_t P11ex_;
      nest::double_t P21ex_;
      nest::double_t P11in_;
      nest::double_t P21in_;
      nest::double_t P22_;
      nest::double_t P30_;
      nest::int_t    RefractoryCounts_;
      nest::int_t    ImmediateInhibitionCounts_;

      librandom::NormalRandomDev normal_dev_;  //!< random deviate generator
    };

    // Access functions for UniversalDataLogger -------------------------------

    //! Read out the real membrane potential
    nest::double_t get_V_m_() const { return S_.V_m_ + P_.E_L_; }

    // ---------------------------------------------------------------- 

    /**
     * @defgroup iaf_psc_exp_data
     * Instances of private data structures for the different types
     * of data pertaining to the model.
     * @note The order of definitions is important for speed.
     * @{
     */   
    Parameters_ P_;
    Buffers_    B_;
    State_      S_;
    Variables_  V_;
    /** @} */

    //! Mapping of recordables names to access functions
    static nest::RecordablesMap<spike_response_neuron> recordablesMap_;  
    
  };

  inline
    nest::port spike_response_neuron::send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex, bool)
  {
    nest::SpikeEvent e;
    e.set_sender(*this);

    return target.handles_test_event(e, receptor_type);
  }

  inline
  nest::port spike_response_neuron::handles_test_event(nest::SpikeEvent&, nest::port receptor_type)
  {
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
 
  inline
  nest::port spike_response_neuron::handles_test_event(nest::DataLoggingRequest &dlr,
					nest::port receptor_type)
  {
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return B_.logger_.connect_logging_device(dlr, recordablesMap_); 
  }

  inline
  void spike_response_neuron::get_status(DictionaryDatum &d) const
  {
    baseparam_.get(d);
    P_.get(d);
    S_.get(d);
    my_base_neuron::get_status(d);

    (*d)[nest::names::recordables] = recordablesMap_.get_list();
  }

  inline
  void spike_response_neuron::set_status(const DictionaryDatum &d)
  {
    BaseParam_ bptmp = baseparam_;            // temporary copy in case of errors
    bptmp.set(d);         // throws if BadProperty
    Parameters_ ptmp = P_;  // temporary copy in case of errors
    ptmp.set(d, S_);
    State_      stmp = S_;  // temporary copy in case of errors
    stmp.set(d);                 // throws if BadProperty

    // We now know that (ptmp, stmp) are consistent. We do not 
    // write them back to (P_, S_) before we are also sure that 
    // the properties to be set in the parent class are internally 
    // consistent.
    my_base_neuron::set_status(d);

    // if we get here, temporaries contain consistent set of properties
    baseparam_ = bptmp;
    P_ = ptmp;
    S_ = stmp;
  }

} // namespace

#endif // SPIKE_RESPONSE_NEURON_H
