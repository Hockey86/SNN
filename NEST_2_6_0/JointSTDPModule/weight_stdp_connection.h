/*
 *  weight_stdp_connection.h
 *
 *  This file is part of JointSTDPModule.
 *
 */

#ifndef WEIGHT_STDP_CONNECTION_H
#define WEIGHT_STDP_CONNECTION_H

/* BeginDocumentation
  Name: weight_stdp_synapse -  Plastic synapse which updates its weight
  according to spike-timing dependent plasticity [1]. It models axonal delays
  and homeostatic nonhebbian terms.

  Description:
    weight_stdp_synapse addictively updates its weight at postsynaptic site.
    In contrast, stdp_synpase in NEST implementation models dendritic delay
    and updates weight at presynaptic site.

  Parameters:
   tau_plus>0    double - Time constant of weight STDP window (positive amp. side), in ms 
   tau_minus>0   double - Time constant of weight STDP window (negative amp. side), in ms 
   A_plus>=0      double - Amplitude of of weight STDP window (positive amp. side)
   A_minus>=0     double - Amplitude of of weight STDP window (negative amp. side)
   w_in>=0       double - Amount of weight change when receiving a presynaptic spike
   w_out<=0      double - Amount of weight change when firing a spike
   STDP_weight_effective_len>0 double - Effective weight STDP window length STDP_weight_effective_len*tau_plus/minus


  Transmits: SpikeEvent
   
  References:
   [3] Song, S., Miller, K. D. and Abbott, L. F. (2000). Competitive
       Hebbian learning through spike-timing-dependent synaptic
       plasticity,Nature Neuroscience 3:9,919--926

  Author: Moritz Helias, Abigail Morrison

  SeeAlso: joint_stdp_synapse, my_base_connection

  Modified: Nov 2015, Sun, Haoqi
*/

#include <cmath>

//#include "connection.h" // included in my_base_connection.h

//#include "network.h"
//#include "dictdatum.h"
//#include "connector_model.h"
//#include "common_synapse_properties.h"
//#include "event.h"
//#include "my_base_neuron.h"
#include "my_base_connection.h"


namespace mynest
{
  template<typename targetidentifierT>
  class JointSTDPConnection;

  // connections are templates of target identifier type (used for pointer / target index addressing)
  // derived from generic connection template
  template<typename targetidentifierT>
  class WeightSTDPConnection : public MyBaseConnection<targetidentifierT>
  {

  public:

    typedef nest::CommonSynapseProperties CommonPropertiesType;
    typedef MyBaseConnection<targetidentifierT> ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
    WeightSTDPConnection();


  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
    WeightSTDPConnection(const WeightSTDPConnection &);

  // Explicitly declare all methods inherited from the dependent base ConnectionBase.
  // This avoids explicit name prefixes in all places these functions are used.
  // Since ConnectionBase depends on the template parameter, they are not automatically
  // found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_delay;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;
  using ConnectionBase::update_weight;
  using ConnectionBase::DT_PRECISION;
  using ConnectionBase::weight_;
  using ConnectionBase::Wmax_;  
  using ConnectionBase::nearest_approximation_;
  using ConnectionBase::debug_;

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status(DictionaryDatum & d) const;

  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status(const DictionaryDatum & d, nest::ConnectorModel &cm);

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param t_lastspike Point in time of last spike sent.
   * \param cp common properties of all synapses (empty).
   */
  void send(nest::Event& e, nest::thread t, nest::double_t t_lastspike, const CommonPropertiesType &cp);
  
  void update_source_upstream(my_base_neuron*, nest::double_t, std::size_t);
  
  void update_delay(nest::double_t) {}
 
  private:

  //nest::double_t facilitate_(nest::double_t w, nest::double_t kplus)
  //{
  //  //return w + lambda_ * std::pow(1.0 - (w/Wmax_), mu_plus_) * kplus * Wmax_;
  //  return w + lambda_ * kplus * Wmax_;
  //}

  //nest::double_t depress_(nest::double_t w, nest::double_t kminus)
  //{
  //  //return w - alpha_ * lambda_ * std::pow(w/Wmax_, mu_minus_) * kminus * Wmax_;
  //  return w - alpha_ * lambda_ * kminus * Wmax_;
  //}
  
  // data members of each connection
  nest::double_t tau_plus_;
  nest::double_t tau_minus_;
  nest::double_t A_plus_;
  nest::double_t A_minus_;
  //nest::double_t mu_plus_;
  //nest::double_t mu_minus_;  
  //nest::double_t Kplus_;
  nest::double_t w_in_;
  nest::double_t w_out_;
  nest::double_t STDP_weight_effective_len_;
  };


  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param t The thread on which this connection is stored.
   * \param t_lastspike Time point of last spike emitted
   * \param cp Common properties object, containing the stdp parameters.
   */
  template<typename targetidentifierT>
  inline
  void WeightSTDPConnection<targetidentifierT>::send(nest::Event& e, nest::thread t, nest::double_t t_lastspike, const CommonPropertiesType &)
  {
    /* downstream stuff */
    
    // get source and target
    my_base_neuron* mytarget = static_cast<my_base_neuron*>(get_target(t));
    my_base_neuron* mysource = static_cast<my_base_neuron*>(&(e.get_sender()));
    if (mytarget==NULL || mysource==NULL)
    {
      std::cout<<std::endl<<std::endl<<"Connection 'weight_stdp_synapse' can only be used together with neurons inherited from my_base_neuron!"<<std::endl<<std::endl;
      throw nest::BadProperty("Connection 'weight_stdp_synapse' can only be used together with neurons inherited from my_base_neuron!");
    }
    
    std::deque<nest::histentry>::iterator start;
    std::deque<nest::histentry>::iterator finish;
    nest::double_t t_spike = e.get_stamp().get_ms(); // t_lastspike_ = 0 initially
    nest::double_t dt, tmp1 = get_delay(), tmp2 = -tau_plus_*STDP_weight_effective_len_;

    std::size_t converted_sender_id = mytarget->convert_upstream_id(mysource->get_gid()), cc;
    /*if (debug_) {
      std::cout << "downstream" << std::endl;
      std::cout << "from neuron #" << mysource->get_gid() << " to #" << mytarget->get_gid() << " original weight: " << weight_ << std::endl;
    }*/
    
    // tell downstream neurons my spike history (I'm the upstream of downstream neurons)
    if (mytarget->get_upstream_connection(converted_sender_id) == NULL) {
      //if (debug_) std::cout<<"target #"<<mytarget->get_gid()<<"does not have upstream connection, create."<<std::endl;
      mytarget->set_upstream_connection(converted_sender_id, this);
    }
    mytarget->append_upstream_spike_history(converted_sender_id, t_spike, t_spike+tmp2-tmp1);
    /*if (debug_) {
      std::list<nest::double_t> ttt = mytarget->get_upstream_spike_history(converted_sender_id);
      for (std::list<nest::double_t>::iterator it = ttt.begin() ; it != ttt.end(); ++it)
        std::cout<<*it<<std::endl;
    }*/
    		      
    // update downstream weights, dW
    // use mytarget->get_history to get downstream neuron spike history
    mytarget->get_history(t_spike+tmp1-STDP_weight_effective_len_*tau_minus_, t_spike, &start, &finish);
    while (start != finish)
    {
      --finish;
      dt = t_spike - finish->t_ + tmp1;
      dt = dt < 0 ? std::ceil(dt*DT_PRECISION-0.5)/DT_PRECISION : std::floor(dt*DT_PRECISION+0.5)/DT_PRECISION;
      //if (debug_) std::cout << t_spike << "-" << finish->t_ << "+" << tmp1 << "=" << dt << std::endl;
      //if (dt > 0.0 && dt <= STDP_weight_effective_len_*tau_minus_)  // dt must be within (delay, STDP_weight_effective_len_*tau_minus_]
        weight_ -= (A_minus_*std::exp(-dt/tau_minus_));
      //else  // dt cannot be within [-STDP_weight_effective_len_*tau_plus_, 0], so comment out
      //  weight_ += (A_plus_*std::exp(dt/tau_plus_));
      if (nearest_approximation_) break;
    }
    //if (debug_) std::cout << "after downstream dW: " << weight_ << std::endl;

    // first order non-Hebbian term on weight, w_in
    update_weight(w_in_);
    //if (debug_) std::cout << "after downstream w_in and cut: " << weight_ << std::endl;
    
    
    /* upstream stuff */
    
    std::size_t source_upstream_num = mysource->get_upstream_num();
    if (source_upstream_num>0)
    {
      //if (debug_) std::cout << "upstream" << std::endl;
      
      // mark this spike is updated, so that different downstream connctions will not update its upstream multiple times
      nest::delay spike_step = e.get_stamp().get_steps();
      if (!mysource->find_updated_spike(spike_step))  // not updated yet
      {        
        mysource->mark_updated_spike(spike_step, spike_step-3);  // mark it as updated, keep a short memory of 3 steps
      
        ConnectionBase* myconn;
          for (cc = 0; cc < source_upstream_num; ++cc)
        {
          if (mysource->get_upstream_connection_type(cc)==0)  // weight_stdp_synapse
            myconn = static_cast<WeightSTDPConnection<targetidentifierT>*>(mysource->get_upstream_connection(cc));
          else  // joint_stdp_synapse
            myconn = static_cast<JointSTDPConnection<targetidentifierT>*>(mysource->get_upstream_connection(cc));
          if (myconn!=NULL)
            myconn->update_source_upstream(mysource, t_spike, cc);
        }
      }
    }

    // send event
    e.set_receiver(*mytarget);
    e.set_weight(weight_);
    e.set_delay(get_delay_steps());
    if (get_delay_steps() <= 0)
    {
      std::cout<<std::endl;
      std::cout<<"=========from weight stdp connection=========: "<<std::endl;
      std::cout<<"delay: "<<e.get_delay()<<std::endl;
      std::cout<<"delay steps: "<<get_delay_steps()<<std::endl;
      std::cout<<"weight: "<<e.get_weight()<<std::endl;
      std::cout<<"receiver gid: "<<e.get_receiver().get_gid()<<std::endl;
      std::cout<<"sender gid: "<<e.get_sender_gid()<<std::endl;
      std::cout<<"stamp [ms] (spike_time): "<<e.get_stamp().get_ms()<<std::endl;
      std::cout<<"max delay: "<<e.get_max_delay()<<std::endl;
      
      std::cout<<"get_upstream_neuron_ids: "<<std::endl;
      std::size_t un = mytarget->get_upstream_num();
      for (std::size_t ii=0;ii<un;++ii)
          std::cout<<mytarget->get_upstream_neuron_id(ii)<<std::endl;
    }
    e.set_rport(get_rport());
    e();
    //if (debug_) std::cout << std::endl;
  }

  template<typename targetidentifierT>
  void WeightSTDPConnection<targetidentifierT>::update_source_upstream(my_base_neuron* source, nest::double_t t_spike, std::size_t cc)
  {
    std::list<nest::double_t> upstream_spike_history = source->get_upstream_spike_history(cc);
    nest::double_t dt, tmp1 = w_out_, tmp2 = tau_minus_*STDP_weight_effective_len_;
    nest::double_t tmp3 = get_delay(), tmp4 = -tau_plus_*STDP_weight_effective_len_;
    
    for (std::list<nest::double_t>::reverse_iterator rit = upstream_spike_history.rbegin() ; rit != upstream_spike_history.rend(); ++rit)
    {
      dt = *rit - t_spike + tmp3;
      dt = dt < 0 ? std::ceil(dt*DT_PRECISION-0.5)/DT_PRECISION : std::floor(dt*DT_PRECISION+0.5)/DT_PRECISION;
      // if outside STDP effective length, break
      // correct only if spike history is sorted ascendly
      //if (dt < tmp4) break;
      // spikes outside STDP effective length is already removed in mytarget->append_upstream_spike_history
      if (dt >= 0.0 && dt <= tmp2)
        tmp1 -= (A_minus_*std::exp(-dt/tau_minus_));
      else if(dt < 0.0 && dt >= tmp4)
        tmp1 += (A_plus_*std::exp(dt/tau_plus_));
      if (debug_) std::cout << cc << ": " << *rit << "-" << t_spike << "+" << tmp3 << "=" << dt << ", dW=" << tmp1-w_out_ << std::endl;
      if (nearest_approximation_) break;
    }
    update_weight(tmp1);
  }


  template<typename targetidentifierT>
  WeightSTDPConnection<targetidentifierT>::WeightSTDPConnection() :
    ConnectionBase(),
    tau_plus_(20.0),
    tau_minus_(20.0),
    A_plus_(0.01),
    A_minus_(0.02),
    //mu_plus_(1.0),
    //mu_minus_(1.0),
    //Kplus_(0.0),
    w_in_(0.0),
    w_out_(0.0),
    STDP_weight_effective_len_(3.0)
  { }

  template<typename targetidentifierT>
  WeightSTDPConnection<targetidentifierT>::WeightSTDPConnection(const WeightSTDPConnection<targetidentifierT> &rhs) :
    ConnectionBase(rhs),
    tau_plus_(rhs.tau_plus_),
    tau_minus_(rhs.tau_minus_),
    A_plus_(rhs.A_plus_),
    A_minus_(rhs.A_minus_),
    //mu_plus_(rhs.mu_plus_),
    //mu_minus_(rhs.mu_minus_),
    //Kplus_(rhs.Kplus_),
    w_in_(rhs.w_in_),
    w_out_(rhs.w_out_),
    STDP_weight_effective_len_(rhs.STDP_weight_effective_len_)
  {  }

  template<typename targetidentifierT>
  void WeightSTDPConnection<targetidentifierT>::get_status(DictionaryDatum & d) const
  {
    ConnectionBase::get_status(d);
    def<nest::double_t>(d, "tau_plus", tau_plus_);
    def<nest::double_t>(d, "tau_minus", tau_minus_);
    def<nest::double_t>(d, "A_plus", A_plus_);
    def<nest::double_t>(d, "A_minus", A_minus_);
    //def<nest::double_t>(d, "mu_plus", mu_plus_);
    //def<nest::double_t>(d, "mu_minus", mu_minus_);
    def<nest::double_t>(d, "w_in", w_in_);
    def<nest::double_t>(d, "w_out", w_out_);
    def<nest::double_t>(d, "STDP_weight_effective_len", STDP_weight_effective_len_);
    def<nest::long_t>(d, nest::names::size_of, sizeof(*this));
  }

  template<typename targetidentifierT>
  void WeightSTDPConnection<targetidentifierT>::set_status(const DictionaryDatum & d, nest::ConnectorModel &cm)
  {
    ConnectionBase::set_status(d, cm);
    updateValue<nest::double_t>(d, "tau_plus", tau_plus_);
    updateValue<nest::double_t>(d, "tau_minus", tau_minus_);
    updateValue<nest::double_t>(d, "A_plus", A_plus_);
    updateValue<nest::double_t>(d, "A_minus", A_minus_);
    //updateValue<nest::double_t>(d, "mu_plus", mu_plus_);
    //updateValue<nest::double_t>(d, "mu_minus", mu_minus_);
    updateValue<nest::double_t>(d, "w_in", w_in_);
    updateValue<nest::double_t>(d, "w_out", w_out_);
    updateValue<nest::double_t>(d, "STDP_weight_effective_len", STDP_weight_effective_len_);
    
    if (tau_plus_ <= 0.0) throw nest::BadProperty("tau_plus should be > 0.");
    if (tau_minus_ <= 0.0) throw nest::BadProperty("tau_minus should be > 0.");
    if (A_plus_ < 0.0) throw nest::BadProperty("A_plus should be >= 0.");
    if (A_minus_ < 0.0) throw nest::BadProperty("A_minus should be >= 0.");
    if (w_out_ > 0.0) throw nest::BadProperty("w_out should be <= 0.");
    if (w_in_ < 0.0) throw nest::BadProperty("w_in should be >= 0.");
    if (STDP_weight_effective_len_ <= 0.0) throw nest::BadProperty("STDP_weight_effective_len should be > 0.");
  }

} // of namespace mynest

#endif // weight_stdp_connection.h
