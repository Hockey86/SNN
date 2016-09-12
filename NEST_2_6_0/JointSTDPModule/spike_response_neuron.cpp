/*
 *  spike_response_neuron.cpp
 *
 *  This file is part of JointSTDPModule.
 *
 */

#include "exceptions.h"
#include "spike_response_neuron.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include "universal_data_logger_impl.h"

#include <limits>

/* ---------------------------------------------------------------- 
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap<mynest::spike_response_neuron> mynest::spike_response_neuron::recordablesMap_;

namespace nest
{
  // Override the create() method with one call to RecordablesMap::insert_() 
  // for each quantity to be recorded.
  template <>
  void RecordablesMap<mynest::spike_response_neuron>::create()
  {
    // use standard names wherever you can for consistency!
    insert_(names::V_m, &mynest::spike_response_neuron::get_V_m_);
  }
}

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
namespace mynest
{

spike_response_neuron::Parameters_::Parameters_()
  : tau_epsp_ (  	8.5        ),  // in ms
    tau_ipsp_ (  	8.5        ),  // in ms
    tau_reset_ (  15.4        ),  // in ms
    E_L_ (   0.0        ),  // normalized
    U_th_ (   1.0				 ),  // normalized
    U_epsp_  (  	0.77       ),  // normalized
    U_ipsp_  (  	0.77       ),  // normalized
    U_reset_  (  	2.31			 ),  // normalized
	C_		 (   1.0				 ),  // Should not be modified
	U_noise_ (   0.0				 ),  // normalized
	noise_	 (						),
	immediate_inhibition_duration_ ( 0.0 ),
	flush_PSP_after_firing_ ( false ),
    TauR_  (  2.0    )  // ms
{}


spike_response_neuron::State_::State_()
  : i_syn_ex_ ( 0.0 ),
    i_syn_in_ ( 0.0 ),
    V_syn_ ( 0.0 ),
    V_spike_ ( 0.0 ),
    V_input_ (0.0),
	V_m_ ( 0.0 ),
    r_  (0)
{}

/* ---------------------------------------------------------------- 
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void spike_response_neuron::Parameters_::get(DictionaryDatum &d) const
{
  def<nest::double_t>(d, nest::names::V_reset,    U_reset_);
  def<nest::double_t>(d, nest::names::V_epsp,     U_epsp_);
  def<nest::double_t>(d, "V_ipsp",     U_ipsp_);
  def<nest::double_t>(d, nest::names::tau_epsp,   tau_epsp_);
  def<nest::double_t>(d, "tau_ipsp",   tau_ipsp_);
  def<nest::double_t>(d, nest::names::tau_reset,  tau_reset_);
  def<nest::double_t>(d, nest::names::V_noise,    U_noise_);
  def<nest::double_t>(d, nest::names::t_ref, TauR_);
  def<nest::double_t>(d, "immediate_inhibition_duration", immediate_inhibition_duration_);
  def<bool>(d, "flush_PSP_after_firing", flush_PSP_after_firing_);
  (*d)[nest::names::noise] = DoubleVectorDatum(new std::vector<double>(noise_));

}

void spike_response_neuron::Parameters_::set(const DictionaryDatum &d, State_& s)
{
  updateValue<nest::double_t>(d, nest::names::V_reset,		U_reset_);
  updateValue<nest::double_t>(d, nest::names::V_epsp, 		U_epsp_);
  updateValue<nest::double_t>(d, "V_ipsp", 		U_ipsp_);
  updateValue<nest::double_t>(d, nest::names::tau_epsp,   tau_epsp_);
  updateValue<nest::double_t>(d, "tau_ipsp",   tau_ipsp_);
  updateValue<nest::double_t>(d, nest::names::tau_reset,  tau_reset_);
  updateValue<nest::double_t>(d, nest::names::V_noise,    U_noise_);
  updateValue<nest::double_t>(d, nest::names::t_ref, TauR_);
  updateValue<nest::double_t>(d, nest::names::V_noise,    U_noise_);
  updateValue<bool>(d, "flush_PSP_after_firing", flush_PSP_after_firing_);
  updateValue<nest::double_t>(d, "immediate_inhibition_duration", immediate_inhibition_duration_);

  const bool updated_noise = updateValue<std::vector<nest::double_t> >(d, nest::names::noise, noise_);
  if(updated_noise) {
  	s.position_ = 0;
  }
  /*
  // TODO: How to handle setting U_noise first and noise later and still make sure they are consistent?
  if ( U_noise_ > 0 && noise_.empty() )
  	throw BadProperty("Noise amplitude larger than zero while noise signal is missing.");
	*/
  if ( U_epsp_ < 0)
    throw nest::BadProperty("V_epsp cannot be negative.");
    
  if ( U_ipsp_ < 0)
    throw nest::BadProperty("V_ipsp cannot be negative. Please use its absolute value.");

  if ( U_reset_ < 0)
    throw nest::BadProperty("Reset potential cannot be negative."); //sign switched above
      
  if ( tau_epsp_ <= 0 ||tau_ipsp_ <= 0 || tau_reset_ <= 0 )
    throw nest::BadProperty("All time constants must be strictly positive.");
}

void spike_response_neuron::State_::get(DictionaryDatum &d) const
{
  def<double>(d, nest::names::V_m, V_m_); // Membrane potential
}

void spike_response_neuron::State_::set(DictionaryDatum const &d)
{
	updateValue<double>(d, nest::names::V_m, V_m_);
}

spike_response_neuron::Buffers_::Buffers_(spike_response_neuron &n)
  : logger_(n)
{}

spike_response_neuron::Buffers_::Buffers_(const Buffers_ &, spike_response_neuron &n)
  : logger_(n)
{}

/* ---------------------------------------------------------------- 
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

spike_response_neuron::spike_response_neuron()
  : my_base_neuron(), 
    P_(), 
    S_(),
    B_(*this)
{
  recordablesMap_.create();
}

spike_response_neuron::spike_response_neuron(const spike_response_neuron &n)
  : my_base_neuron(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void spike_response_neuron::init_node_(const nest::Node &proto)
{
  const spike_response_neuron &pr = downcast<spike_response_neuron>(proto);
  P_ = pr.P_;
  S_ = pr.S_;
}

void spike_response_neuron::init_state_(const nest::Node &proto)
{
  const spike_response_neuron &pr = downcast<spike_response_neuron>(proto);
  S_ = pr.S_;
}

void spike_response_neuron::init_buffers_()
{
  B_.spikes_ex_.clear();        // includes resize
  B_.spikes_in_.clear();        // includes resize
  B_.currents_.clear();         // includes resize
  B_.logger_.reset();
  my_base_neuron::clear_history();
  my_base_neuron::clear_upstream_spike_history();
}

void spike_response_neuron::calibrate()
{
  B_.logger_.init();  // ensures initialization in case mm connected after Simulate

  const double h = nest::Time::get_resolution().get_ms(); 

  // numbering of state vaiables: i_0 = 0, i_syn_ = 1, V_syn_ = 2, V_spike _= 3, V_m_ = 4

  // these P are independent
  V_.P11ex_ = std::exp(-h/P_.tau_epsp_);
  V_.P11in_ = std::exp(-h/P_.tau_ipsp_);

  V_.P22_ = std::exp(-h/P_.tau_epsp_);

  V_.P30_ = std::exp(-h/P_.tau_reset_);

  // these depend on the above. Please do not change the order.
  // TODO: use expm1 here to improve accuracy for small timesteps

  V_.P21ex_ = P_.U_epsp_*std::exp(1.0)/(P_.C_) * V_.P11ex_ * h/P_.tau_epsp_;
  V_.P21in_ = P_.U_ipsp_*std::exp(1.0)/(P_.C_) * V_.P11in_ * h/P_.tau_ipsp_;

  V_.P20_ = P_.tau_epsp_/P_.C_*(1.0 - V_.P22_);
  
  V_.RefractoryCounts_ = nest::Time(nest::Time::ms(P_.TauR_)).get_steps();
  V_.ImmediateInhibitionCounts_ = nest::Time(nest::Time::ms(P_.immediate_inhibition_duration_)).get_steps();
  assert(V_.RefractoryCounts_ >= 0);  // since t_ref_ >= 0, this can only fail in error
  assert(V_.ImmediateInhibitionCounts_ >= 0);
}

void spike_response_neuron::update(const nest::Time &origin, const nest::long_t from, const nest::long_t to)
{
  assert(to >= 0 && (nest::delay) from < nest::Scheduler::get_min_delay());
  assert(from < to);

  // evolve from timestep 'from' to timestep 'to' with steps of h each
  for ( nest::long_t lag = from; lag < to; ++lag )
  {	
    S_.V_syn_ = S_.V_syn_*V_.P22_ + S_.i_syn_ex_*V_.P21ex_ + S_.i_syn_in_*V_.P21in_;
      
    // exponential decaying PSCs
    S_.i_syn_ex_ *= V_.P11ex_;
    S_.i_syn_in_ *= V_.P11in_;

    // the spikes arriving at T+1 have an immediate effect on the state of the neuron
    S_.i_syn_ex_ += B_.spikes_ex_.get_value(lag);
    S_.i_syn_in_ += B_.spikes_in_.get_value(lag);

    // exponentially decaying ahp
    S_.V_spike_ *= V_.P30_;
    
    double noise_term = P_.U_noise_ > 0.0 && !P_.noise_.empty() ? P_.U_noise_ * P_.noise_[S_.position_++] : 0.0;

    if ( S_.r_ == 0 )
    {
      // neuron not refractory
      S_.V_m_ = S_.V_syn_ + S_.V_spike_ + noise_term + S_.V_input_;
    }
    else // neuron is absolute refractory
      --S_.r_;
      
    // cleaned immediately after use
    if (S_.V_input_>0.0) S_.V_input_ = 0.0;

                                                       
    if ( S_.V_m_ >= P_.U_th_ )  // threshold crossing
    {
      S_.r_  = V_.RefractoryCounts_;
      //S_.V_spike_ -= P_.U_reset_;
      //S_.V_m_ -= P_.U_reset_;
      S_.V_spike_ = -P_.U_reset_;
      S_.V_m_ = -P_.U_reset_;
      
      if ( P_.flush_PSP_after_firing_)
      {
        S_.i_syn_ex_ = 0.0;
        S_.i_syn_in_ = 0.0;
        S_.V_syn_ = 0.0;
      }

      set_spiketime(nest::Time::step(origin.get_steps()+lag+1));
	    
      nest::SpikeEvent se;
      network()->send(*this, se, lag);
      
      for (int i=0; i<baseparam_.downstream_immediate_inhibition_neuron_num_; ++i)
      {
        if (baseparam_.downstream_immediate_inhibition_neurons_[i]!=NULL)
        {
          //nest::SpikeEvent* myse = new nest::SpikeEvent(se);
          //myse->set_weight(baseparam_.downstream_weights_[it-baseparam_.downstream_immediate_inhibition_neurons_.begin()]);
          //myse->set_delay(0);
          //(*it)->handle(*myse);
          //delete myse;            
          baseparam_.downstream_immediate_inhibition_neurons_[i]->immediate_inhibit(baseparam_.downstream_weights_[i]);
        }
      }
    }

    // log state data
    B_.logger_.record_data(origin.get_steps() + lag);

  }  
}                           
                     
void spike_response_neuron::handle(nest::SpikeEvent &e)
{
  if (e.get_delay() <= 0)
  {
    std::cout<<std::endl;
    std::cout<<"=========from spike_response_neuron=========: "<<std::endl;
    std::cout<<"delay: "<<e.get_delay()<<std::endl;
    std::cout<<"weight: "<<e.get_weight()<<std::endl;
    std::cout<<"multiplicity: "<<e.get_multiplicity()<<std::endl;
    std::cout<<"rel_delivery_steps: "<<e.get_rel_delivery_steps(network()->get_slice_origin())<<std::endl;
    std::cout<<"receiver gid: "<<e.get_receiver().get_gid()<<std::endl;
    std::cout<<"sender gid: "<<e.get_sender_gid()<<std::endl;
    std::cout<<"stamp [ms]: "<<e.get_stamp().get_ms()<<std::endl;
    std::cout<<"max delay: "<<e.get_max_delay()<<std::endl;
  }
  assert ( e.get_delay() > 0 );

  if ( e.get_weight() > 0.0 )
    B_.spikes_ex_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
			    e.get_weight() * e.get_multiplicity() );
  else
    B_.spikes_in_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
           e.get_weight() * e.get_multiplicity() );
}

void spike_response_neuron::handle(nest::DataLoggingRequest &e)
{
  B_.logger_.handle(e);
}
} // namespace
