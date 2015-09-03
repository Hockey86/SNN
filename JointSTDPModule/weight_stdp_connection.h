/*
 *  weight_stdp_connection.h
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

#ifndef WEIGHT_STDP_CONNECTION_H
#define WEIGHT_STDP_CONNECTION_H

/* BeginDocumentation
  Name: weight_stdp_synapse - Synapse type for (additive) spike-timing
  dependent plasticity with non-Hebbian first-order terms.

  Description:
   weight_stdp_synapse is a connector to create synapses with spike time 
   dependent plasticity (as defined in [1]) with non-Hebbian first-order terms.
   Here the weight dependence exponent can be set separately
   for potentiation and depression.

  Examples:
   //multiplicative STDP [2]  mu_plus = mu_minus = 1.0
   additive STDP       [3]  mu_plus = mu_minus = 0.0
   //Guetig STDP         [1]  mu_plus = mu_minus = [0.0,1.0]
   //van Rossum STDP     [4]  mu_plus = 0.0 mu_minus = 1.0 

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


  Transmits: SpikeEvent
   
  References:
   [1] Guetig et al. (2003) Learning Input Correlations through Nonlinear
       Temporally Asymmetric Hebbian Plasticity. Journal of Neuroscience

   [2] Rubin, J., Lee, D. and Sompolinsky, H. (2001). Equilibrium
       properties of temporally asymmetric Hebbian plasticity, PRL
       86,364-367

   [3] Song, S., Miller, K. D. and Abbott, L. F. (2000). Competitive
       Hebbian learning through spike-timing-dependent synaptic
       plasticity,Nature Neuroscience 3:9,919--926

   [4] van Rossum, M. C. W., Bi, G-Q and Turrigiano, G. G. (2000). 
       Stable Hebbian learning from spike timing-dependent
       plasticity, Journal of Neuroscience, 20:23,8812--8821

  FirstVersion: Fri 00:37 May 22 2015 @ NTU
  Author: Haoqi Sun

  MODIFIED FROM
  Author: Moritz Helias, Abigail Morrison

  SeeAlso: synapsedict, tsodyks_synapse, static_synapse
*/

#include <cmath>

#include "connection.h"

//#include "network.h"
//#include "dictdatum.h"
//#include "connector_model.h"
//#include "common_synapse_properties.h"
//#include "event.h"
#include "iaf_psc_alpha2.h"


namespace mynest
{

  // connections are templates of target identifier type (used for pointer / target index addressing)
  // derived from generic connection template
  template<typename targetidentifierT>
  class WeightSTDPConnection : public nest::Connection<targetidentifierT>
  {

  public:

    typedef nest::CommonSynapseProperties CommonPropertiesType;
    typedef nest::Connection<targetidentifierT> ConnectionBase;

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

  void set_weight(nest::double_t w) { weight_ = w; }
  
  void update_weight(nest::double_t dw) { weight_ += dw;
  if (weight_ >= Wmax_) weight_ = Wmax_;
  else if (weight_ <= 0.0) weight_ = 0.0;
  }
 
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
void WeightSTDPConnection<targetidentifierT>::send(nest::Event& e, nest::thread t, nest::double_t t_lastspike, const CommonPropertiesType &)
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
    throw nest::BadProperty("Connection type 'weight_stdp_synapse' can only be used together with neuron type 'iaf_psc_alpha2' !");
  }
  
  std::deque<nest::histentry>::iterator start;
  std::deque<nest::histentry>::iterator finish;
  std::vector<nest::double_t> spike_history, upstream_spike_history;
  nest::double_t t_spike = e.get_stamp().get_ms(); // t_lastspike_ = 0 initially
  nest::double_t dt, tmp1 = get_delay(), tmp2, tmp3, tmp4;

  nest::index sender_gid_1 = mysource->get_gid()-1;
  nest::index target_gid_1 = mytarget->get_gid()-1;
  std::size_t converted_sender_id = mytarget->convert_id(sender_gid_1), cc;
  //if (debug_) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " original weight: " << weight_ << std::endl;
  
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
		      
  // update downstream weights, dW
  //tmp2 = STDP_weight_window_plus_*tau_plus_;
  mytarget->get_history(t_spike+tmp1-STDP_weight_window_plus_*tau_plus_, t_spike, &start, &finish);
  //std::cout << "downstream" << std::endl;
  while (start != finish)
  {
    //if (debug_) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " downstream spike time: " << start->t_ << std::endl;
    dt = t_spike - start->t_ + tmp1;
    //std::cout << dt << std::endl;
    //if (dt > 0.0 && dt <= tmp2)  // dt must be within (delay, STDP_weight_window_plus_*tau_plus_]
    //{
      weight_ -= (A_plus_*std::exp(-dt/tau_plus_));
      /*if (STDP_window_type_==0) weight_ -= (A_plus_*std::exp(-dt/tau_plus_));
      else if (STDP_window_type_==1) weight_ -= (2.71828182846*A_plus_/tau_plus_*dt*std::exp(-dt/tau_plus_));
      else if (STDP_window_type_==2) weight_ -= (1.648721*A_plus_/tau_plus_*dt*std::exp(-std::pow(dt/tau_plus_,2)/2.0));*/
    //}
    //else  // dt cannot be within [-STDP_weight_window_minus_*tau_minus_, 0], so comment out
    //{
    //  weight_ += (A_minus_*std::exp(dt/tau_minus_));
      /*if (STDP_window_type_==0) weight_ += (A_minus_*std::exp(dt/tau_minus_));
      else if (STDP_window_type_==1) weight_ -= (2.71828182846*A_minus_/tau_minus_*dt*std::exp(dt/tau_minus_));
      else if (STDP_window_type_==2) weight_ -= (1.648721*A_minus_/tau_minus_*dt*std::exp(-std::pow(dt/tau_minus_,2)/2.0));*/
    //}
    ++start;
  }
  //if (debug_) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " after downstream dW: " << weight_ << std::endl;

  // first order non-Hebbian term on weight, a_pre
  weight_ += a_pre_;
  //if (debug_) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " after downstream a_pre: " << weight_ << std::endl;
  
  // cut the weight
  if (weight_ >= Wmax_) weight_ = Wmax_;
  else if (weight_ <= 0.0) weight_ = 0.0;
  //if (debug_) std::cout << "from neuron #" << sender_gid_1 << " to #" << target_gid_1 << " after cut: " << weight_ << std::endl;
  
  
  
  /* upstream stuff */
  
  std::size_t source_upstream_num = mysource->get_upstream_num();
  if (source_upstream_num>0)
  {
    WeightSTDPConnection* myconn;
    tmp2 = tau_plus_*STDP_weight_window_plus_;
    tmp3 = -tau_minus_*STDP_weight_window_minus_;
    //std::cout << "upstream" << std::endl;
    for (cc = 0; cc < source_upstream_num; ++cc)
    {
      myconn = static_cast<WeightSTDPConnection*>(mysource->get_upstream_connection(cc));
      if (myconn!=NULL)
      {
        tmp1 = a_post_;        
        tmp4 = myconn->get_delay();
        upstream_spike_history = mysource->get_upstream_spike_history(cc);
        for (std::vector<nest::double_t>::iterator it = upstream_spike_history.begin() ; it != upstream_spike_history.end(); ++it)
        {
          dt = *it - t_spike + tmp4;
          //std::cout << cc << " " << dt << std::endl;
          if (dt > 0.0 && dt <= tmp2)
          {
            tmp1 -= (A_plus_*std::exp(-dt/tau_plus_));
            /*if (STDP_window_type_==0) tmp1 -= (A_plus_*std::exp(-dt/tau_plus_));
            else if (STDP_window_type_==1) tmp1 -= (2.71828182846*A_plus_/tau_plus_*dt*std::exp(-dt/tau_plus_));
            else if (STDP_window_type_==2) tmp1 -= (1.648721*A_plus_/tau_plus_*dt*std::exp(-std::pow(dt/tau_plus_,2)/2.0));*/
          }
          else if(dt <= 0.0 && dt >= tmp3)
          {
            tmp1 += (A_minus_*std::exp(dt/tau_minus_));
            /*if (STDP_window_type_==0) tmp1 += (A_minus_*std::exp(dt/tau_minus_));
            else if (STDP_window_type_==1) tmp1 -= (2.71828182846*A_minus_/tau_minus_*dt*std::exp(dt/tau_minus_));
            else if (STDP_window_type_==2) tmp1 -= (1.648721*A_minus_/tau_minus_*dt*std::exp(-std::pow(dt/tau_minus_,2)/2.0));*/
          }
        }
        myconn->update_weight(tmp1);
      }
    }
  }

  // send event
  e.set_receiver(*mytarget);
  e.set_weight(weight_);
  e.set_delay(get_delay_steps());
  e.set_rport(get_rport());
  e();
  //if (debug_) std::cout << std::endl;
}


  template<typename targetidentifierT>
  WeightSTDPConnection<targetidentifierT>::WeightSTDPConnection() :
    ConnectionBase(),
    weight_(1.0),
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
    STDP_window_type_(0),
    debug_(false)
  { }

  template<typename targetidentifierT>
  WeightSTDPConnection<targetidentifierT>::WeightSTDPConnection(const WeightSTDPConnection<targetidentifierT> &rhs) :
    ConnectionBase(rhs),
    weight_(rhs.weight_),
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
    STDP_window_type_(rhs.STDP_window_type_),
    debug_(rhs.debug_)
  {  }

  template<typename targetidentifierT>
  void WeightSTDPConnection<targetidentifierT>::get_status(DictionaryDatum & d) const
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
    def<nest::long_t>(d, "STDP_window_type", STDP_window_type_);
    def<bool>(d, "debug", debug_);
    def<nest::long_t>(d, nest::names::size_of, sizeof(*this));
  }

  template<typename targetidentifierT>
  void WeightSTDPConnection<targetidentifierT>::set_status(const DictionaryDatum & d, nest::ConnectorModel &cm)
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
    updateValue<nest::long_t>(d, "STDP_window_type", STDP_window_type_);
    updateValue<bool>(d, "debug", debug_);
  }

} // of namespace mynest

#endif // weight_stdp_connection.h
