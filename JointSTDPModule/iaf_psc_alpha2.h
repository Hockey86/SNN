/*
 *  iaf_psc_alpha2.h
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

#ifndef IAF_PSC_ALPHA2_H
#define IAF_PSC_ALPHA2_H

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"
#include "recordables_map.h"

#include "nestmodule.h"
#include "target_identifier.h"

/* BeginDocumentation
Name: iaf_psc_alpha2 - Leaky integrate-and-fire neuron model, with storage of upstream delays.

Description:

  iaf_psc_alpha2 is an implementation of a leaky integrate-and-fire model
  with alpha-function shaped synaptic currents. Thus, synaptic currents
  and the resulting post-synaptic potentials have a finite rise time. 

  The threshold crossing is followed by an absolute refractory period
  during which the membrane potential is clamped to the resting potential.

  The linear subthresold dynamics is integrated by the Exact
  Integration scheme [1]. The neuron dynamics is solved on the time
  grid given by the computation step size. Incoming as well as emitted
  spikes are forced to that grid.  

  An additional state variable and the corresponding differential
  equation represents a piecewise constant external current.

  The general framework for the consistent formulation of systems with
  neuron like dynamics interacting by point events is described in
  [1].  A flow chart can be found in [2].

  Critical tests for the formulation of the neuron model are the
  comparisons of simulation results for different computation step
  sizes. sli/testsuite/nest contains a number of such tests.
  
  The iaf_psc_alpha2 is the standard model used to check the consistency
  of the nest simulation kernel because it is at the same time complex
  enough to exhibit non-trivial dynamics and simple enough compute
  relevant measures analytically.

Remarks:

  The present implementation uses individual variables for the
  components of the state vector and the non-zero matrix elements of
  the propagator.  Because the propagator is a lower triangular matrix
  no full matrix multiplication needs to be carried out and the
  computation can be done "in place" i.e. no temporary state vector
  object is required.

  The template support of recent C++ compilers enables a more succinct
  formulation without loss of runtime performance already at minimal
  optimization levels. A future version of iaf_psc_alpha2 will probably
  address the problem of efficient usage of appropriate vector and
  matrix objects.


Parameters: 

  The following parameters can be set in the status dictionary.

  V_m        double - Membrane potential in mV 
  E_L        double - Resting membrane potential in mV. 
  C_m        double - Capacity of the membrane in pF
  tau_m      double - Membrane time constant in ms.
  t_ref      double - Duration of refractory period in ms. 
  V_th       double - Spike threshold in mV.
  V_reset    double - Reset potential of the membrane in mV.
  tau_syn_ex double - Rise time of the excitatory synaptic alpha function in ms.
  tau_syn_in double - Rise time of the inhibitory synaptic alpha function in ms.
  I_e        double - Constant external input current in pA.
  V_min      double - Absolute lower value for the membrane potential.
 
Note:
  tau_m != tau_syn_{ex,in} is required by the current implementation to avoid a
  degenerate case of the ODE describing the model [1]. For very similar values,
  numerics will be unstable.

References:
  [1] Rotter S & Diesmann M (1999) Exact simulation of time-invariant linear
      systems with applications to neuronal modeling. Biologial Cybernetics
      81:381-402.
  [2] Diesmann M, Gewaltig M-O, Rotter S, & Aertsen A (2001) State space 
      analysis of synchronous spiking in cortical neural networks. 
      Neurocomputing 38-40:565-571.
  [3] Morrison A, Straube S, Plesser H E, & Diesmann M (2006) Exact subthreshold 
      integration with continuous spike times in discrete time neural network 
      simulations. Neural Computation, in press

Sends: nest::SpikeEvent

Receives: nest::SpikeEvent, nest::CurrentEvent, nest::DataLoggingRequest

FirstVersion: Mon 22:53 May 25 2015 @ NTU
Author:  Haoqi Sun

MODIFIED FROM
FirstVersion: September 1999
Author:  Diesmann, Gewaltig
SeeAlso: iaf_psc_delta, iaf_psc_exp, iaf_cond_exp
*/

namespace mynest
{
  class Network;

  /**
   * Leaky integrate-and-fire neuron with alpha-shaped PSCs.
   */
  class iaf_psc_alpha2 : public nest::Archiving_Node
  {
    
  public:
    
    iaf_psc_alpha2();
    iaf_psc_alpha2(const iaf_psc_alpha2&);

    /**
     * Import sets of overloaded virtual functions.
     * @see Technical Issues / Virtual Functions: Overriding, Overloading, and Hiding
     */
    using nest::Node::handle;
    using nest::Node::handles_test_event;
    bool has_proxies()    const { return false; }
    bool local_receiver() const { return true;  }

    nest::port send_test_event(nest::Node&, nest::rport, nest::synindex, bool);
    
    void handle(nest::SpikeEvent &);
    void handle(nest::CurrentEvent &);
    void handle(nest::DataLoggingRequest &); 
    
    nest::port handles_test_event(nest::SpikeEvent&, nest::rport);
    nest::port handles_test_event(nest::CurrentEvent&, nest::rport);
    nest::port handles_test_event(nest::DataLoggingRequest &, nest::rport);

    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);

    std::vector<nest::Connection<nest::TargetIdentifierPtrRport>* > get_upstream_connections()
    { return P_.upstream_connections_; }

    nest::Connection<nest::TargetIdentifierPtrRport>* get_upstream_connection(nest::size_t id)
    { return P_.upstream_connections_[id]; }

    std::vector<nest::long_t> get_upstream_neuron_ids() const
    { return P_.upstream_neuron_ids_; }

    //std::vector<nest::double_t> get_upstream_spike_times() const
    //{ return P_.upstream_spike_times_; }

    std::vector<nest::double_t> get_upstream_delays() const
    { return P_.upstream_delays_; }

    //nest::double_t get_upstream_spike_time(std::size_t id) const
    //{ return P_.upstream_spike_times_[id]; }

    std::vector<nest::double_t> get_upstream_spike_history(std::size_t id) const
    { return P_.upstream_spike_history_[id]; }

    nest::double_t get_upstream_delay(std::size_t id) const
    { return P_.upstream_delays_[id]; }

    std::size_t convert_id(nest::size_t id)
    {
      std::vector<nest::long_t>::iterator loc = std::find( P_.upstream_neuron_ids_.begin(), P_.upstream_neuron_ids_.end(), id);
      if ( loc != P_.upstream_neuron_ids_.end() )
        return (std::size_t)(loc-P_.upstream_neuron_ids_.begin());
      else
        throw nest::BadProperty("Upstream id not found.");
    }
    
    std::size_t get_upstream_num() { return P_.upstream_neuron_ids_.size(); }

    void set_upstream_connection(nest::size_t id, nest::Connection<nest::TargetIdentifierPtrRport>* conn)
    { P_.upstream_connections_[id] = conn; }

    void set_upstream_spike_history(nest::size_t id, std::vector<nest::double_t> ush)
    { P_.upstream_spike_history_[id] = ush; }
    
    //void set_upstream_spike_time(nest::size_t id, nest::double_t tf)
    //{ P_.upstream_spike_times_[id] = tf; }

    void set_upstream_delay(nest::size_t id, nest::double_t dv)
    { P_.upstream_delays_[id] = dv; }

  private:

    void init_state_(const nest::Node& proto);
    void init_buffers_();
    void calibrate();

    void update(nest::Time const &, const nest::long_t, const nest::long_t);

    // The next two classes need to be friends to access the State_ class/member
    friend class nest::RecordablesMap<iaf_psc_alpha2>;
    friend class nest::UniversalDataLogger<iaf_psc_alpha2>;

    // ---------------------------------------------------------------- 

    struct Parameters_ {

      /** upstream connections */
      std::vector<nest::Connection<nest::TargetIdentifierPtrRport>*> upstream_connections_;
      
      /** upstream delays */
      std::vector<nest::double_t> upstream_delays_;

      /** upstream spike history */
      //std::vector<nest::double_t> upstream_spike_times_;
      std::vector<std::vector<nest::double_t> > upstream_spike_history_;

      /** upstream neuron ids */
      std::vector<nest::long_t> upstream_neuron_ids_;
  
      /** Membrane time constant in ms. */
      nest::double_t Tau_; 

      /** Membrane capacitance in pF. */
      nest::double_t C_;
    
      /** Refractory period in ms. */
      nest::double_t TauR_;

      /** Resting potential in mV. */
      nest::double_t U0_;

      /** External current in pA */
      nest::double_t I_e_;

      /** Reset value of the membrane potential */
      nest::double_t V_reset_;

      /** Threshold, RELATIVE TO RESTING POTENTIAL(!).
          I.e. the real threshold is (U0_+Theta_). */
      nest::double_t Theta_;

      /** Lower bound, RELATIVE TO RESTING POTENTIAL(!).
          I.e. the real lower bound is (LowerBound_+U0_). */
      nest::double_t LowerBound_;

      /** nest::Time constant of excitatory synaptic current in ms. */
      nest::double_t tau_ex_;

      /** nest::Time constant of inhibitory synaptic current in ms. */
      nest::double_t tau_in_;
      
      Parameters_();  //!< Sets default parameter values

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary

      /** Set values from dictionary.
       * @returns Change in reversal potential E_L, to be passed to State_::set()
       */
      double set(const DictionaryDatum&);

    };
    
    // ---------------------------------------------------------------- 

    struct State_ {

      nest::double_t y0_; //!< Constant current
      nest::double_t y1_ex_;  
      nest::double_t y2_ex_;
      nest::double_t y1_in_;
      nest::double_t y2_in_;
      nest::double_t y3_; //!< This is the membrane potential RELATIVE TO RESTING POTENTIAL.

      nest::int_t    r_;  //!< Number of refractory steps remaining

      State_();  //!< Default initialization
      
      void get(DictionaryDatum&, const Parameters_&) const;

      /** Set values from dictionary.
       * @param dictionary to take data from
       * @param current parameters
       * @param Change in reversal potential E_L specified by this dict
       */
      void set(const DictionaryDatum&, const Parameters_&, double);

    };

    // ---------------------------------------------------------------- 

    struct Buffers_ {

      Buffers_(iaf_psc_alpha2&);
      Buffers_(const Buffers_&, iaf_psc_alpha2&);

      /** buffers and summs up incoming spikes/currents */
      nest::RingBuffer ex_spikes_;
      nest::RingBuffer in_spikes_;
      nest::RingBuffer currents_;

      //! Logger for all analog data
      nest::UniversalDataLogger<iaf_psc_alpha2> logger_;

    };
    
    // ---------------------------------------------------------------- 

    struct Variables_ {

      /** Amplitude of the synaptic current.
	  This value is chosen such that a post-synaptic potential with
	  weight one has an amplitude of 1 mV.
       */
      nest::double_t EPSCInitialValue_;
      nest::double_t IPSCInitialValue_;
      nest::int_t    RefractoryCounts_;
    
      nest::double_t P11_ex_;
      nest::double_t P21_ex_;
      nest::double_t P22_ex_;
      nest::double_t P31_ex_;
      nest::double_t P32_ex_;
      nest::double_t P11_in_;
      nest::double_t P21_in_;
      nest::double_t P22_in_;
      nest::double_t P31_in_;
      nest::double_t P32_in_;
      nest::double_t P30_;
      nest::double_t P33_;
      nest::double_t expm1_tau_m_;

      nest::double_t weighted_spikes_ex_;
      nest::double_t weighted_spikes_in_;


    };

    // Access functions for UniversalDataLogger -------------------------------

    //! Read out the real membrane potential
    nest::double_t get_V_m_() const { return S_.y3_ + P_.U0_; }

    nest::double_t get_weighted_spikes_ex_() const { return V_.weighted_spikes_ex_; }
    nest::double_t get_weighted_spikes_in_() const { return V_.weighted_spikes_in_; }
    nest::double_t get_input_currents_ex_() const { return S_.y1_ex_; }
    nest::double_t get_input_currents_in_() const { return S_.y1_in_; }

    // Data members ----------------------------------------------------------- 
    
    /**
     * @defgroup iaf_psc_alpha2_data
     * Instances of private data structures for the different types
     * of data pertaining to the model.
     * @note The order of definitions is important for speed.
     * @{
     */   
    Parameters_ P_;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;
    /** @} */
    
    //! Mapping of recordables nest::names to access functions
    static nest::RecordablesMap<iaf_psc_alpha2> recordablesMap_;
  };

  inline
  nest::port iaf_psc_alpha2::send_test_event(nest::Node& target, nest::rport receptor_type, nest::synindex, bool)
  {
    nest::SpikeEvent e;
    e.set_sender(*this);
    return target.handles_test_event(e, receptor_type);
  }
    
  inline
  nest::port iaf_psc_alpha2::handles_test_event(nest::SpikeEvent&, nest::rport receptor_type)
  {
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
   
  inline
  nest::port iaf_psc_alpha2::handles_test_event(nest::CurrentEvent&, nest::rport receptor_type)
  {
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
   
  inline
  nest::port iaf_psc_alpha2::handles_test_event(nest::DataLoggingRequest& dlr, 
  				       nest::rport receptor_type)
  {
    if (receptor_type != 0)
      throw nest::UnknownReceptorType(receptor_type, get_name());
    return B_.logger_.connect_logging_device(dlr, recordablesMap_);
  }

  inline
  void iaf_psc_alpha2::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    S_.get(d, P_);
    nest::Archiving_Node::get_status(d);
  
    (*d)[nest::names::recordables] = recordablesMap_.get_list();
  }
  
  inline
  void iaf_psc_alpha2::set_status(const DictionaryDatum &d)
  {
    Parameters_ ptmp = P_;            // temporary copy in case of errors
    const double delta_EL = ptmp.set(d);         // throws if BadProperty
    State_      stmp = S_;            // temporary copy in case of errors
    stmp.set(d, ptmp, delta_EL);                 // throws if BadProperty
  
    // We now know that (ptmp, stmp) are consistent. We do not 
    // write them back to (P_, S_) before we are also sure that 
    // the properties to be set in the parent class are internally 
    // consistent.
    nest::Archiving_Node::set_status(d);
  
    // if we get here, temporaries contain consistent set of properties
    P_ = ptmp;
    S_ = stmp;
  }

} // namespace

#endif /* #ifndef IAF_PSC_ALPHA2_H */
