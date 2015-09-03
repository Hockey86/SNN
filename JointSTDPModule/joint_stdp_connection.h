/*
 *  joint_stdp_connection.h
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

#ifndef JOINT_STDP_CONNECTION_H
#define JOINT_STDP_CONNECTION_H

/* BeginDocumentation
  Name: joint_stdp_synapse - Synapse type for additive joint
  weight-delay spike-timing dependent plasticity..

  Description:
   joint_stdp_synapse is a connector to create synapses with spike time 
   dependent plasticity (as defined in [1]) on both weight and delay.

  Parameters:
   tau_plus    double - Time constant of weight STDP window (dt=pre-post>0 side), in ms 
   tau_minus   double - Time constant of weight STDP window (dt=pre-post<=0 side), in ms 
   A_plus      double - Amplitude of of weight STDP window (dt=pre-post>0 side)
   A_minus     double - Amplitude of of weight STDP window (dt=pre-post<=0 side)
   Wmax        double - Maximum allowed weight
   a_pre       double - Amount of weight change when receiving a presynaptic spike
   a_post      double - Amount of weight change when firing a spike
   STDP_weight_window_plus double - Valid window length for weight update, STDP_weight_window_plus*tau_plus
   STDP_weight_window_minus  double - Valid window length for weight update, STDP_weight_window_minus*tau_minus
   
   sigma_plus  double - Time constant of delay STDP window (dt=pre-post>0 side), in ms 
   sigma_minus double - Time constant of delay STDP window (dt=pre-post<=0 side), in ms 
   B_plus      double - Amplitude of of delay STDP window (dt=pre-post>0 side)
   B_minus     double - Amplitude of of delay STDP window (dt=pre-post<=0 side)
   STDP_delay_window_plus double - Valid window length for delay update, STDP_delay_window_plus*sigma_plus
   STDP_delay_window_minus  double - Valid window length for delay update, STDP_delay_window_minus*sigma_minus
   delay_min   double - minimum delay in ms


  Transmits: SpikeEvent
   
  References:
   [1] Guetig et al. (2003) Learning Input Correlations through Nonlinear
       Temporally Asymmetric Hebbian Plasticity. Journal of Neuroscience

  FirstVersion: Fri 00:37 May 22 2015 @ NTU
  Author: Haoqi Sun

  MODIFIED FROM
  Author: Moritz Helias, Abigail Morrison

  SeeAlso: synapsedict, tsodyks_synapse, static_synapse
*/

#include <cmath>
#include <limits>

#include "connection.h"

//#include "network.h"
//#include "dictdatum.h"
//#include "connector_model.h"
//#include "common_synapse_properties.h"
//#include "event.h"
#include "nest_time.h"
#include "iaf_psc_alpha2.h"


namespace mynest
{
  // connections are templates of target identifier type (used for pointer / target index addressing)
  // derived from generic connection template
  template<typename targetidentifierT>
  class JointSTDPConnection : public nest::Connection<targetidentifierT>
  {

  public:

    typedef nest::CommonSynapseProperties CommonPropertiesType;
    typedef nest::Connection<targetidentifierT> ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
    JointSTDPConnection();


  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
    JointSTDPConnection(const JointSTDPConnection &);

  // Explicitly declare all methods inherited from the dependent base ConnectionBase.
  // This avoids explicit name prefixes in all places these functions are used.
  // Since ConnectionBase depends on the template parameter, they are not automatically
  // found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_delay;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;
  using ConnectionBase::set_delay_steps;

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


  class ConnTestDummyNode: public nest::ConnTestDummyNodeBase
  {
  public:
	// Ensure proper overriding of overloaded virtual functions.
	// Return values from functions are ignored.
	using nest::ConnTestDummyNodeBase::handles_test_event;
    nest::port handles_test_event(nest::SpikeEvent&, nest::rport) { return nest::invalid_port_; }
  };

  void check_connection(nest::Node & s, nest::Node & t, nest::rport receptor_type, nest::double_t t_lastspike, const CommonPropertiesType &)
  {
    ConnTestDummyNode dummy_target;
    
    ConnectionBase::check_connection_(dummy_target, s, t, receptor_type);

    t.register_stdp_connection(t_lastspike - get_delay());
  }

  nest::double_t get_weight() { return weight_; }

  void set_weight(nest::double_t w) { weight_ = w; }
  
  void update_weight(nest::double_t dw) { weight_ += dw;
  if (weight_ > Wmax_) weight_ = Wmax_;
  else if (weight_ < 0.0) weight_ = 0.0;
  }
  
 protected:
  // ConnectionBase already provides 'delay', but is rounded to time_step
  // we need to store the precise delay to avoid rounding error
  
  // *internal_delay_* equals axonal delay + forward dendritic delay
  // for simplicity we assume axonal delay = forward dendritic delay = backpropagating dendritic delay
  // we have *internal_delay_* = 2*axonal delay
  // and the different between pre- and post-synaptic firing time is t_i^f-t^f+ (axonal delay-backpropagating dendritic delay) = t_i^f-t^f
  nest::double_t internal_delay_;
 
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
  nest::double_t weight_;
  nest::double_t Wmax_;
  nest::double_t tau_plus_;
  nest::double_t tau_minus_;
  nest::double_t A_plus_;
  nest::double_t A_minus_;
  //nest::double_t mu_plus_;
  //nest::double_t mu_minus_;  
  //nest::double_t Kplus_;
  nest::double_t a_pre_;
  nest::double_t a_post_;
  nest::double_t STDP_weight_window_plus_;
  nest::double_t STDP_weight_window_minus_;
  
  //nest::double_t eta_d_;
  nest::double_t delay_min_;
  nest::double_t sigma_plus_;
  nest::double_t sigma_minus_;
  nest::double_t B_plus_;
  nest::double_t B_minus_;
  nest::double_t STDP_delay_window_plus_;
  nest::double_t STDP_delay_window_minus_;
  nest::long_t STDP_window_type_;
  bool debug_;

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
void JointSTDPConnection<targetidentifierT>::send(nest::Event& e, nest::thread t, nest::double_t t_lastspike, const CommonPropertiesType &)
{
  /* downstream stuff */
  
  // get source and target
  iaf_psc_alpha2* mytarget = static_cast<iaf_psc_alpha2*>(get_target(t));
  iaf_psc_alpha2* mysource = static_cast<iaf_psc_alpha2*>(&(e.get_sender()));
  if (mytarget==NULL || mysource == NULL)
  {
    //mytarget = static_cast<iaf_neuron2*>(target);
    //mysource = static_cast<iaf_neuron2&>(e.get_sender());
    //if (mytarget==NULL || mysource==NULL)
    throw nest::BadProperty("Connection type 'joint_stdp_synapse' can only be used together with neuron type 'iaf_psc_alpha2' !");
  }
  
  std::deque<nest::histentry>::iterator start;
  std::deque<nest::histentry>::iterator finish;
  std::vector<nest::double_t> spike_history, upstream_spike_history;
  nest::double_t t_spike = e.get_stamp().get_ms(); // t_lastspike_ = 0 initially
  nest::double_t dt, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6;
  if (internal_delay_ < 0.0)
    internal_delay_ = get_delay();
  nest::double_t delay_max_ = nest::Time::delay_steps_to_ms(e.get_max_delay());
  //nest::double_t original_dendritic_delay = dendritic_delay;

  nest::index sender_gid_1 = mysource->get_gid()-1;
  nest::index target_gid_1 = mytarget->get_gid()-1;
  std::size_t converted_sender_id = mytarget->convert_id(sender_gid_1), cc;
  //if (debug_) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " original weight: " << weight_ << " original delay: " << internal_delay_ << std::endl;
  
  // tell downstream neurons my spike history (I'm the upstream of downstream neurons)
  mysource->get_history(t_spike - STDP_weight_window_minus_*tau_minus_, t_spike-0.01, &start, &finish);
  while (start != finish)
  {
    spike_history.push_back(start->t_);
    ++start;
  }
  spike_history.push_back(t_spike);
  mytarget->set_upstream_spike_history(converted_sender_id, spike_history);
  if (mytarget->get_upstream_connection(converted_sender_id) == NULL)
    mytarget->set_upstream_connection(converted_sender_id, this);
		      
  // update downstream weights and delays
  //tmp1 = t_spike-STDP_delay_window_minus_*sigma_minus_;
  tmp2 = STDP_weight_window_plus_*tau_plus_;
  mytarget->get_history(t_spike+internal_delay_-STDP_weight_window_plus_*tau_plus_, t_spike, &start, &finish);
  //std::cout << "downstream" << std::endl;
  while (start != finish)
  {
    //if (debug_) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " downstream spike time: " << start->t_ << std::endl;
    dt = t_spike - start->t_ + internal_delay_;
    //std::cout << dt << std::endl;
    //if (dt > 0.0 && dt <= tmp2)  // dt must be within (delay, STDP_weight_window_plus_*tau_plus_]
    //{
    weight_ -= (A_plus_*std::exp(-dt/tau_plus_));
      /*if (STDP_window_type_==0) weight_ -= (A_plus_*std::exp(-dt/tau_plus_));
      else if (STDP_window_type_==1) weight_ -= (2.71828182846*A_plus_/tau_plus_*dt*std::exp(-dt/tau_plus_));
      else if (STDP_window_type_==2) weight_ -= (1.648721*A_plus_/tau_plus_*dt*std::exp(-std::pow(dt/tau_plus_,2)/2.0));*/
      //if (debug_) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " weight changed to: " << weight_ << " because I fired at " << t_spike <<" downstream neuron fired at " << start->t_ << std::endl;
    //}
    //else  // dt cannot be within [-STDP_weight_window_minus_*tau_minus_, 0], so comment out
    //{
    //  weight_ += (A_minus_*std::exp(dt/tau_minus_));
      /*if (STDP_window_type_==0) weight_ += (A_minus_*std::exp(dt/tau_minus_));
      else if (STDP_window_type_==1) weight_ -= (2.71828182846*A_minus_/tau_minus_*dt*std::exp(dt/tau_minus_));
      else if (STDP_window_type_==2) weight_ -= (1.648721*A_minus_/tau_minus_*dt*std::exp(-std::pow(dt/tau_minus_,2)/2.0));*/
      //if (debug_) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " weight changed to: " << weight_ << " because I fired at " << t_spike <<" downstream neuron fired at " << start->t_ << std::endl;
    //}
    ++start;
  }
  tmp4 = STDP_delay_window_plus_*sigma_plus_;
  mytarget->get_history(t_spike+internal_delay_-STDP_weight_window_plus_*sigma_plus_, t_spike, &start, &finish);
  while (start != finish)
  {
    dt = t_spike - start->t_ + internal_delay_;
    //if (dt > 0.0 && dt <= tmp4)  // dt must be within (delay, STDP_weight_window_plus_*sigma_plus_]
    //{
    internal_delay_ += (B_plus_*std::exp(-dt/sigma_plus_));
      /*if (STDP_window_type_==0) internal_delay_ += (B_plus_*std::exp(-dt/sigma_plus_));
      else if (STDP_window_type_==1) internal_delay_ += (2.71828182846*B_plus_/sigma_plus_*dt*std::exp(-dt/sigma_plus_));
      else if (STDP_window_type_==2) internal_delay_ += (1.648721*B_plus_/sigma_plus_*dt*std::exp(-std::pow(dt/sigma_plus_,2)/2.0));*/
      //if (debug_) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " delay changed to: " << internal_delay_ << " because I fired at " << t_spike <<" downstream neuron fired at " << start->t_ << std::endl;
    //}
    //else  // dt cannot be within [-STDP_weight_window_minus_*sigma_minus_, 0], so comment out
    //  internal_delay_ += ?? (B_minus_*std::exp(dt/sigma_minus_));
    ++start;
  }
  //if (debug_) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " after downstream dW: " << weight_ << std::endl;

  // first order non-Hebbian term on weight, a_pre
  weight_ += a_pre_;
  //if (debug_) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " after downstream a_pre: " << weight_ << std::endl;
  
  // cut the weight
  //if (debug_ && weight_ > Wmax_) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " weight cut to: " << Wmax_ << std::endl;
  //if (debug_ && weight_ < 0.0) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " weight cut to: " << 0.0 << std::endl;
  if (weight_ > Wmax_) weight_ = Wmax_;
  else if (weight_ < 0.0) weight_ = 0.0;
  //if (debug_) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " after cut: " << weight_ << std::endl;
  
  // cut the delay
  //if (debug_ && internal_delay_ > delay_max_) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " delay cut to: " << delay_max_ << std::endl;
  //if (debug_ && internal_delay_ <  delay_min_) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " delay cut to: " <<  delay_min_ << std::endl;
  if (internal_delay_ > delay_max_) internal_delay_ = delay_max_;
  else if (internal_delay_ < delay_min_) internal_delay_ = delay_min_;
  
  
  
  /* upstream stuff */
  
  std::size_t source_upstream_num = mysource->get_upstream_num();
  if (source_upstream_num>0)
  {
  
    // Delta d_ij = -( min_k (d_ik+tkf) - (d_ij+tjf) )* SIGMA_q {B_minus_*exp((upstream_spike_q-t_spike)/sigma_minus_)}
    // first, get upstream delays from upstream connections (if NULL, from user input)
    std::vector<ConnectionBase*> upstream_connections = mysource->get_upstream_connections();
    std::vector<nest::double_t> upstream_spike_times, upstream_delays;
    //std::vector<nest::long_t> upstream_neuron_ids;///////
    //if (debug_) std::cout << "neuron #" << sender_gid_1 << " upstream delays: "<< std::endl;
    for(typename std::vector<ConnectionBase*>::iterator it = upstream_connections.begin() ; it != upstream_connections.end(); ++it)
    {
      if ((*it) == NULL)
        upstream_delays.push_back(mysource->get_upstream_delay(it-upstream_connections.begin()));
      else
        upstream_delays.push_back((static_cast<JointSTDPConnection*>(*it))->internal_delay_);
      //if (debug_) std::cout << upstream_delays.back() << std::endl;
    }
    // second, get upstream spike times
    for (cc = 0; cc < source_upstream_num; cc++)
    {
      if (mysource->get_upstream_spike_history(cc).size()>0)
        upstream_spike_times.push_back(mysource->get_upstream_spike_history(cc).back());
      else
        upstream_spike_times.push_back(-1.0); // no spike, spike time is negative
    }
    
    tmp3 = t_spike-sigma_minus_*STDP_delay_window_minus_;

    // third, get earliest presynaptic arrival time min_k (d_ik+tkf)
    nest::double_t earliest_presynaptic_arrival = std::numeric_limits<nest::double_t>::max();
    std::vector<bool> has_presynaptic_spike;
    //if (debug_) upstream_neuron_ids = mysource->get_upstream_neuron_ids();
    //if (debug_) std::cout << "neuron #" << sender_gid_1 << " upstream spike times: "<< std::endl;
    for (std::vector<nest::double_t>::iterator it = upstream_spike_times.begin() ; it != upstream_spike_times.end(); ++it)
    {
      cc = it-upstream_spike_times.begin();
      //if (debug_) std::cout << "#" << upstream_neuron_ids[cc] << ": ";
      
      tmp1 = *it+upstream_delays[cc];
      
      if (*it>0.0 && *it >= tmp3 && tmp1 <= t_spike) // future event tmp1 > t_spike is not allowed (but deserve better solution)
      {
        has_presynaptic_spike.push_back(true);
        //if (debug_) std::cout << *it << " + " << upstream_delays[cc] << " = " << tmp1 << std::endl;
        if (tmp1 < earliest_presynaptic_arrival)
          earliest_presynaptic_arrival = tmp1;
        //if (debug_) std::cout << *it << std::endl;
      }
      else
      {
        has_presynaptic_spike.push_back(false);
        /*if (debug_)
        {
          if (*it<=0.0)
            std::cout << "no presynaptic spike" << std::endl;
          else if (tmp1 > t_spike)
            std::cout << "this spike arrives later than postsynatic firing: " << tmp1 << ", now: " << t_spike << std::endl;
          else
            std::cout << "out of window, presynaptic spike: " << *it << ", now: " << t_spike << std::endl;
        }*/
      }
    }
    //if (debug_) std::cout << "neuron #" << sender_gid_1 << " earliest spike arrival: "<< earliest_presynaptic_arrival << std::endl;
    
    // now we know the earliest presynaptic arrival time min_k (d_ik+tkf)
    // make upstream changes on delay and weight
    JointSTDPConnection* myconn;
    tmp4 = sigma_plus_*STDP_delay_window_plus_;
    tmp5 = tau_plus_*STDP_weight_window_plus_;
    tmp6 = -tau_minus_*STDP_weight_window_minus_;
    //std::cout << "upstream" << std::endl;
    for (std::vector<nest::double_t>::iterator it = upstream_delays.begin() ; it != upstream_delays.end(); ++it)
    {
      cc = it-upstream_delays.begin();
      myconn = static_cast<JointSTDPConnection*>(mysource->get_upstream_connection(cc));
      if (myconn!=NULL)
      {
        tmp1 = a_post_;
        tmp2 = 0.0;
        tmp3 = earliest_presynaptic_arrival-*it-upstream_spike_times[cc];
        
        upstream_spike_history = mysource->get_upstream_spike_history(cc);
        for (std::vector<nest::double_t>::iterator itt = upstream_spike_history.begin() ; itt != upstream_spike_history.end(); ++itt)
        {
          dt = *itt - t_spike + *it;
          //std::cout << cc << " " << dt << std::endl;
          if (dt > 0.0 && dt <= tmp5)
          {
            tmp1 -= (A_plus_*std::exp(-dt/tau_plus_));
            /*if (STDP_window_type_==0) tmp1 -= (A_plus_*std::exp(-dt/tau_plus_));
            else if (STDP_window_type_==1) tmp1 -= (2.71828182846*A_plus_/tau_plus_*dt*std::exp(-dt/tau_plus_));
            else if (STDP_window_type_==2) tmp1 -= (1.648721*A_plus_/tau_plus_*dt*std::exp(-std::pow(dt/tau_plus_,2)/2.0));*/
          }
          else if (dt <= 0.0 && dt >= tmp6)
          {
            tmp1 += (A_minus_*std::exp(dt/tau_minus_));
            /*if (STDP_window_type_==0) tmp1 += (A_minus_*std::exp(dt/tau_minus_));
            else if (STDP_window_type_==1) tmp1 -= (2.71828182846*A_minus_/tau_minus_*dt*std::exp(dt/tau_minus_));
            else if (STDP_window_type_==2) tmp1 -= (1.648721*A_minus_/tau_minus_*dt*std::exp(-std::pow(dt/tau_minus_,2)/2.0));*/
          }
          if (dt > 0.0 && dt <= tmp4)
          {
            tmp2 += (B_plus_*std::exp(-dt/sigma_plus_));
            /*if (STDP_window_type_==0) tmp2 += (B_plus_*std::exp(-dt/sigma_plus_));
            else if (STDP_window_type_==1) tmp2 += (2.71828182846*B_plus_/sigma_plus_*dt*std::exp(-dt/sigma_plus_));
            else if (STDP_window_type_==2) tmp2 += (1.648721*B_plus_/sigma_plus_*dt*std::exp(-std::pow(dt/sigma_plus_,2)/2.0));*/
          }
          else if (has_presynaptic_spike[cc])
          {
            tmp2 += (tmp3*B_minus_*std::exp(dt/sigma_minus_));
            /*if (STDP_window_type_==0) tmp2 += (tmp3*B_minus_*std::exp(dt/sigma_minus_));
            else if (STDP_window_type_==1) tmp2 -= (2.71828182846*tmp3*B_minus_/sigma_minus_*dt*std::exp(dt/sigma_minus_));
            else if (STDP_window_type_==2) tmp2 -= (tmp3*1.648721*B_minus_/sigma_minus_*dt*std::exp(-std::pow(dt/sigma_minus_,2)/2.0));*/
          }
        }
        myconn->update_weight(tmp1);
        *it += tmp2;
        //if (debug_) std::cout << "neuron #" << sender_gid_1 << " upstream neuron #" << cc << " weight changed to: " << myconn->get_weight() << " delay changed to: " << *it << " because upstream neuron fired at (see above) " << " I fired at " << t_spike << std::endl;
      }
      /*else
      {      
        if (debug_) std::cout << "neuron #" << sender_gid_1 << " upstream neuron #" << cc << " connection unknown." << std::endl;
      }*/
    }
    //std::cout << std::endl;
    
    /*if (debug_)
    {
      std::cout << "downstream neuron #" << sender_gid_1 << " before shift" << std::endl;
      for (std::vector<nest::double_t>::iterator it = upstream_delays.begin() ; it != upstream_delays.end(); ++it)
        std::cout << "upstream neuron #" << upstream_neuron_ids[it-upstream_delays.begin()] << ": " << *it << std::endl;
    }*/
    
    tmp1 = *std::min_element(upstream_delays.begin(), upstream_delays.end());
    if (tmp1 < delay_min_)
    {
      // shift all upstream delays of the source neurons so that min_i d_i = delay_min_
      for (std::vector<nest::double_t>::iterator it = upstream_delays.begin() ; it != upstream_delays.end(); ++it)
      {
        tmp2 = *it-tmp1+delay_min_;
        *it = tmp2 < delay_max_ ? tmp2 : delay_max_;
      }
    }
    else
    {
      for (std::vector<nest::double_t>::iterator it = upstream_delays.begin() ; it != upstream_delays.end(); ++it)
        if (*it > delay_max_)
          *it = delay_max_;
    }
    /*if (debug_)
    {
      std::cout << "downstream neuron #" << sender_gid_1 << " after shift" << std::endl;
      for (std::vector<nest::double_t>::iterator it = upstream_delays.begin() ; it != upstream_delays.end(); ++it)
        std::cout << "upstream neuron #" << upstream_neuron_ids[it-upstream_delays.begin()] << ": " << *it << std::endl;
    }*/
    
    // update upstream delays
    for (std::vector<nest::double_t>::iterator it = upstream_delays.begin() ; it != upstream_delays.end(); ++it)
    {
      cc = it-upstream_delays.begin();
      myconn = static_cast<JointSTDPConnection*>(mysource->get_upstream_connection(cc));
      if (myconn!=NULL)
      {
        myconn->internal_delay_ = *it;
        myconn->set_delay_steps(nest::Time::delay_ms_to_steps(*it));
        //if (debug_) std::cout << "neuron #" << sender_gid_1 << " upstream neuron #" << cc << " delay cut to: " << *it << std::endl;
      }
    }
  }

  // send event
  e.set_receiver(*mytarget);
  e.set_weight(weight_);
  e.set_delay(nest::Time::delay_ms_to_steps(internal_delay_));
  e.set_rport(get_rport());
  e();
  //if (debug_) std::cout << std::endl;
}


  template<typename targetidentifierT>
  JointSTDPConnection<targetidentifierT>::JointSTDPConnection() :
    ConnectionBase(),
    weight_(1.0),
    internal_delay_(-1.0),
    tau_plus_(20.0),
    tau_minus_(20.0),
    A_plus_(0.01),
    A_minus_(0.02),
    //mu_plus_(1.0),
    //mu_minus_(1.0),
    Wmax_(100.0),
    //Kplus_(0.0),
    a_pre_(0.0),
    a_post_(0.0),
    STDP_weight_window_plus_(3.0),
    STDP_weight_window_minus_(3.0),
    //eta_d_(0.001),
    delay_min_(1.0),
    sigma_plus_(20.0),
    sigma_minus_(20.0),
    B_plus_(0.4),
    B_minus_(0.2),
    STDP_delay_window_plus_(3.0),
    STDP_delay_window_minus_(3.0),
    STDP_window_type_(0),
    debug_(false)
  { }

  template<typename targetidentifierT>
  JointSTDPConnection<targetidentifierT>::JointSTDPConnection(const JointSTDPConnection<targetidentifierT> &rhs) :
    ConnectionBase(rhs),
    weight_(rhs.weight_),
    internal_delay_(rhs.internal_delay_),
    tau_plus_(rhs.tau_plus_),
    tau_minus_(rhs.tau_minus_),
    A_plus_(rhs.A_plus_),
    A_minus_(rhs.A_minus_),
    //mu_plus_(rhs.mu_plus_),
    //mu_minus_(rhs.mu_minus_),
    Wmax_(rhs.Wmax_),
    //Kplus_(rhs.Kplus_),
    a_pre_(rhs.a_pre_),
    a_post_(rhs.a_post_),
    STDP_weight_window_plus_(rhs.STDP_weight_window_plus_),
    STDP_weight_window_minus_(rhs.STDP_weight_window_minus_),
    //eta_d_(rhs.eta_d_),
    delay_min_(rhs.delay_min_),
    sigma_plus_(rhs.sigma_plus_),
    sigma_minus_(rhs.tau_minus_),
    B_plus_(rhs.B_plus_),
    B_minus_(rhs.B_minus_),
    STDP_delay_window_plus_(rhs.STDP_delay_window_plus_),
    STDP_delay_window_minus_(rhs.STDP_delay_window_minus_),
    STDP_window_type_(rhs.STDP_window_type_),
    debug_(rhs.debug_)
  { }

  template<typename targetidentifierT>
  void JointSTDPConnection<targetidentifierT>::get_status(DictionaryDatum & d) const
  {
    ConnectionBase::get_status(d);
    def<nest::double_t>(d, nest::names::weight, weight_);
    def<nest::double_t>(d, "tau_plus", tau_plus_);
    def<nest::double_t>(d, "tau_minus", tau_minus_);
    def<nest::double_t>(d, "A_plus", A_plus_);
    def<nest::double_t>(d, "A_minus", A_minus_);
    //def<nest::double_t>(d, "mu_plus", mu_plus_);
    //def<nest::double_t>(d, "mu_minus", mu_minus_);
    def<nest::double_t>(d, "Wmax", Wmax_);
    def<nest::double_t>(d, "a_pre", a_pre_);
    def<nest::double_t>(d, "a_post", a_post_);
    def<nest::double_t>(d, "STDP_weight_window_plus", STDP_weight_window_plus_);
    def<nest::double_t>(d, "STDP_weight_window_minus", STDP_weight_window_minus_);
    //def<nest::double_t>(d, "eta_d", eta_d_);
    def<nest::double_t>(d, "delay_min", delay_min_);
    def<nest::double_t>(d, "sigma_plus", sigma_plus_);
    def<nest::double_t>(d, "sigma_minus", sigma_minus_);
    def<nest::double_t>(d, "B_plus", B_plus_);
    def<nest::double_t>(d, "B_minus", B_minus_);
    def<nest::double_t>(d, "STDP_delay_window_plus", STDP_delay_window_plus_);
    def<nest::double_t>(d, "STDP_delay_window_minus", STDP_delay_window_minus_);
    def<nest::long_t>(d, "STDP_window_type", STDP_window_type_);
    def<bool>(d, "debug", debug_);
    def<nest::long_t>(d, nest::names::size_of, sizeof(*this));
  }

  template<typename targetidentifierT>
  void JointSTDPConnection<targetidentifierT>::set_status(const DictionaryDatum & d, nest::ConnectorModel &cm)
  {
    ConnectionBase::set_status(d, cm);
    updateValue<nest::double_t>(d, nest::names::weight, weight_);
    updateValue<nest::double_t>(d, "tau_plus", tau_plus_);
    updateValue<nest::double_t>(d, "tau_minus", tau_minus_);
    updateValue<nest::double_t>(d, "A_plus", A_plus_);
    updateValue<nest::double_t>(d, "A_minus", A_minus_);
    //updateValue<nest::double_t>(d, "mu_plus", mu_plus_);
    //updateValue<nest::double_t>(d, "mu_minus", mu_minus_);
    updateValue<nest::double_t>(d, "Wmax", Wmax_);
    updateValue<nest::double_t>(d, "a_pre", a_pre_);
    updateValue<nest::double_t>(d, "a_post", a_post_);
    updateValue<nest::double_t>(d, "STDP_weight_window_plus", STDP_weight_window_plus_);
    updateValue<nest::double_t>(d, "STDP_weight_window_minus", STDP_weight_window_minus_);
    //updateValue<nest::double_t>(d, "eta_d", eta_d_);
    updateValue<nest::double_t>(d, "delay_min", delay_min_);
    updateValue<nest::double_t>(d, "sigma_plus", sigma_plus_);
    updateValue<nest::double_t>(d, "sigma_minus", sigma_minus_);
    updateValue<nest::double_t>(d, "B_plus", B_plus_);
    updateValue<nest::double_t>(d, "B_minus", B_minus_);
    updateValue<nest::double_t>(d, "STDP_delay_window_plus", STDP_delay_window_plus_);
    updateValue<nest::double_t>(d, "STDP_delay_window_minus", STDP_delay_window_minus_);
    updateValue<nest::long_t>(d, "STDP_window_type", STDP_window_type_);
    if (STDP_delay_window_minus_ > STDP_weight_window_minus_)
    {
      std::cout << "STDP_delay_window_minus must be equal or shorter than STDP_weight_window_minus." << std::endl;
      throw nest::BadProperty("STDP_delay_window_minus must be equal or shorter than STDP_weight_window_minus.");
    }
    updateValue<bool>(d, "debug", debug_);
  }

} // of namespace mynest

#endif // joint_stdp_connection.h
