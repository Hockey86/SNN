/*
 *  iaf_psc_alpha2.cpp
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

#include "exceptions.h"
#include "iaf_psc_alpha2.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include "universal_data_logger_impl.h"

#include <limits>

nest::RecordablesMap<mynest::iaf_psc_alpha2> mynest::iaf_psc_alpha2::recordablesMap_;

namespace nest
{

  /*
   * Override the create() method with one call to RecordablesMap::insert_()
   * for each quantity to be recorded.
   */
  template <>
  void RecordablesMap<mynest::iaf_psc_alpha2>::create()
  {
    // use standard nest::names whereever you can for consistency!
    insert_(names::V_m, &mynest::iaf_psc_alpha2::get_V_m_);
    insert_(names::weighted_spikes_ex, &mynest::iaf_psc_alpha2::get_weighted_spikes_ex_);
    insert_(names::weighted_spikes_in, &mynest::iaf_psc_alpha2::get_weighted_spikes_in_);
    insert_(names::input_currents_ex, &mynest::iaf_psc_alpha2::get_input_currents_ex_);
    insert_(names::input_currents_in, &mynest::iaf_psc_alpha2::get_input_currents_in_);
  }
}

namespace mynest
{
  /* ----------------------------------------------------------------
   * Default constructors defining default parameters and state
   * ---------------------------------------------------------------- */

  iaf_psc_alpha2::Parameters_::Parameters_()
    : Tau_       ( 10.0    ),  // ms
      C_         (250.0    ),  // pF
      TauR_      (  2.0    ),  // ms
      U0_        (-70.0    ),  // mV
      I_e_       (  0.0    ),  // pA
      V_reset_   (-70.0-U0_),  // mV, rel to U0_
      Theta_     (-55.0-U0_),  // mV, rel to U0_
      LowerBound_(-std::numeric_limits<nest::double_t>::infinity()),
      tau_ex_    (  2.0    ),  // ms
      tau_in_    (  2.0    )   // ms
  {}

  iaf_psc_alpha2::State_::State_()
    : y0_   (0.0),
      y1_ex_(0.0),
      y2_ex_(0.0),
      y1_in_(0.0),
      y2_in_(0.0),
      y3_   (0.0),
      r_    (0)
  {}

  /* ----------------------------------------------------------------
   * Parameter and state extractions and manipulation functions
   * ---------------------------------------------------------------- */

  void iaf_psc_alpha2::Parameters_::get(DictionaryDatum &d) const
  {
    def<nest::double_t>(d, nest::names::E_L, U0_);   // Resting potential
    def<nest::double_t>(d, nest::names::I_e, I_e_);
    def<nest::double_t>(d, nest::names::V_th, Theta_+U0_); // threshold value
    def<nest::double_t>(d, nest::names::V_reset, V_reset_+U0_);
    def<nest::double_t>(d, nest::names::V_min, LowerBound_+U0_);
    def<nest::double_t>(d, nest::names::C_m, C_);
    def<nest::double_t>(d, nest::names::tau_m, Tau_);
    def<nest::double_t>(d, nest::names::t_ref, TauR_);
    def<nest::double_t>(d, nest::names::tau_syn_ex, tau_ex_);
    def<nest::double_t>(d, nest::names::tau_syn_in, tau_in_);
    def<std::vector<nest::double_t> >(d, "upstream_delays", upstream_delays_);
    def<std::vector<nest::long_t> >(d, "upstream_neuron_ids", upstream_neuron_ids_);
  }

  nest::double_t iaf_psc_alpha2::Parameters_::set(const DictionaryDatum& d)
  {
    // if U0_ is changed, we need to adjust all variables defined relative to U0_
    const nest::double_t ELold = U0_;
    updateValue<nest::double_t>(d, nest::names::E_L, U0_);
    const nest::double_t delta_EL = U0_ - ELold;

    if(updateValue<nest::double_t>(d, nest::names::V_reset, V_reset_))
      V_reset_ -= U0_;
    else
      V_reset_ -= delta_EL;

    if (updateValue<nest::double_t>(d, nest::names::V_th, Theta_))
      Theta_ -= U0_;
    else
      Theta_ -= delta_EL;

    if (updateValue<nest::double_t>(d, nest::names::V_min, LowerBound_))
      LowerBound_ -= U0_;
    else
      LowerBound_ -= delta_EL;

    updateValue<nest::double_t>(d, nest::names::I_e, I_e_);
    updateValue<nest::double_t>(d, nest::names::C_m, C_);
    updateValue<nest::double_t>(d, nest::names::tau_m, Tau_);
    updateValue<nest::double_t>(d, nest::names::tau_syn_ex, tau_ex_);
    updateValue<nest::double_t>(d, nest::names::tau_syn_in, tau_in_);
    updateValue<nest::double_t>(d, nest::names::t_ref, TauR_);
    updateValue<std::vector<nest::double_t> >(d, "upstream_delays", upstream_delays_);
    updateValue<std::vector<nest::long_t> >(d, "upstream_neuron_ids", upstream_neuron_ids_);

    if ( C_ <= 0.0 )
      throw nest::BadProperty("I'm not good at dealing with non-positive capacitance. Ready.. crash!");

    if ( Tau_ <= 0.0 )
      throw nest::BadProperty("Non-positive membrane time constant detected. This type of neuron only exist in the alien brain.");

    if (tau_ex_ <= 0.0 || tau_in_ <= 0.0 )
      throw nest::BadProperty("Neurons with superfast synaptic time constants <= 0 only exist in your penis.");

    if ( Tau_ == tau_ex_ || Tau_ == tau_in_ )
      throw nest::BadProperty("Membrane and synapse time constant(s) must differ. See note in documentation. Don't ask me.");

    if ( TauR_ < 0.0 )
    	throw nest::BadProperty("Just don't know how to deal with negative refractory. Go back to the past?");

    if ( V_reset_ >= Theta_ )
      throw nest::BadProperty("You must be kidding (reset potential higher than threshold) LOL.");

    if ( upstream_delays_.size() != upstream_neuron_ids_.size() )
      throw nest::BadProperty("5 upstream neuron ids have 6 delays? You suck at mathematics.");
      
    if ( upstream_delays_.size() > 0 && *std::min_element(upstream_delays_.begin(), upstream_delays_.end()) < 0.1 )
      throw nest::BadProperty("LOL, non-positive upstream delays!");
      
    if ( upstream_neuron_ids_.size() > 0 && *std::min_element(upstream_neuron_ids_.begin(), upstream_neuron_ids_.end()) < 0 )
      throw nest::BadProperty("I love neuroscience, but I hate negative upstream neurons id.");
      
    upstream_spike_history_ = std::vector<std::vector<nest::double_t> >(upstream_delays_.size(),std::vector<nest::double_t>(0));
    upstream_connections_ = std::vector<nest::Connection<nest::TargetIdentifierPtrRport>* >(upstream_delays_.size(), NULL);

    return delta_EL;
  }

  void iaf_psc_alpha2::State_::get(DictionaryDatum &d, const Parameters_& p) const
  {
    def<nest::double_t>(d, nest::names::V_m, y3_ + p.U0_); // Membrane potential
  }

  void iaf_psc_alpha2::State_::set(const DictionaryDatum& d, const Parameters_& p, nest::double_t delta_EL)
  {
    if ( updateValue<nest::double_t>(d, nest::names::V_m, y3_) )
      y3_ -= p.U0_;
    else
      y3_ -= delta_EL;
  }

  iaf_psc_alpha2::Buffers_::Buffers_(iaf_psc_alpha2& n)
    : logger_(n)
  {}

  iaf_psc_alpha2::Buffers_::Buffers_(const Buffers_ &, iaf_psc_alpha2& n)
    : logger_(n)
  {}


  /* ----------------------------------------------------------------
   * Default and copy constructor for node
   * ---------------------------------------------------------------- */

  iaf_psc_alpha2::iaf_psc_alpha2()
    : nest::Archiving_Node(),
      P_(),
      S_(),
      B_(*this)
  {
    recordablesMap_.create();
  }

  iaf_psc_alpha2::iaf_psc_alpha2(const iaf_psc_alpha2& n)
    : nest::Archiving_Node(n),
      P_(n.P_),
      S_(n.S_),
      B_(n.B_, *this)
  {}

  /* ----------------------------------------------------------------
   * Node initialization functions
   * ---------------------------------------------------------------- */

  void iaf_psc_alpha2::init_state_(const nest::Node& proto)
  {
    const iaf_psc_alpha2& pr = downcast<iaf_psc_alpha2>(proto);
    S_ = pr.S_;
  }

  void iaf_psc_alpha2::init_buffers_()
  {
    B_.ex_spikes_.clear();       // includes resize
    B_.in_spikes_.clear();       // includes resize
    B_.currents_.clear();        // includes resize

    B_.logger_.reset();

    nest::Archiving_Node::clear_history();
    
    for (std::size_t i = 0; i < P_.upstream_spike_history_.size(); i++)
      P_.upstream_spike_history_[i].clear();
  }

  void iaf_psc_alpha2::calibrate()
  {
    B_.logger_.init();  // ensures initialization in case mm connected after Simulate

    const nest::double_t h = nest::Time::get_resolution().get_ms();

    // these P are independent
    V_.P11_ex_ = V_.P22_ex_ = std::exp(-h/P_.tau_ex_);
    V_.P11_in_ = V_.P22_in_ = std::exp(-h/P_.tau_in_);

    V_.P33_ = std::exp(-h/P_.Tau_);

    V_.expm1_tau_m_ = numerics::expm1(-h/P_.Tau_);

    // these depend on the above. Please do not change the order.
    V_.P30_ = -P_.Tau_/P_.C_*numerics::expm1(-h/P_.Tau_);

    V_.P21_ex_ = h * V_.P11_ex_;
    V_.P31_ex_ = 1/P_.C_ * ((V_.P11_ex_-V_.P33_)/(-1/P_.tau_ex_- -1/P_.Tau_)- h*V_.P11_ex_)
      /(-1/P_.Tau_ - -1/P_.tau_ex_);
    V_.P32_ex_ = 1/P_.C_*(V_.P33_-V_.P11_ex_)/(-1/P_.Tau_ - -1/P_.tau_ex_);

    V_.P21_in_ = h * V_.P11_in_;
    V_.P31_in_ = 1/P_.C_ * ((V_.P11_in_-V_.P33_)/(-1/P_.tau_in_- -1/P_.Tau_)- h*V_.P11_in_)
      /(-1/P_.Tau_ - -1/P_.tau_in_);
    V_.P32_in_ = 1/P_.C_*(V_.P33_-V_.P11_in_)/(-1/P_.Tau_ - -1/P_.tau_in_);

    V_.EPSCInitialValue_=1.0 * numerics::e/P_.tau_ex_;
    V_.IPSCInitialValue_=1.0 * numerics::e/P_.tau_in_;

    // TauR specifies the length of the absolute refractory period as
    // a double_t in ms. The grid based iaf_psc_alpha2 can only handle refractory
    // periods that are integer multiples of the computation step size (h).
    // To ensure consistency with the overall simulation scheme such conversion
    // should be carried out via objects of class nest::Time. The conversion
    // requires 2 steps:
    //     1. A time object is constructed defining representation of
    //        TauR in tics. This representation is then converted to computation time
    //        steps again by a strategy defined by class nest::Time.
    //     2. The refractory time in units of steps is read out get_steps(), a member
    //        function of class nest::Time.
    //
    // The definition of the refractory period of the iaf_psc_alpha2 is consistent
    // the one of iaf_psc_alpha2_ps.
    //
    // Choosing a TauR that is not an integer multiple of the computation time
    // step h will lead to accurate (up to the resolution h) and self-consistent
    // results. However, a neuron model capable of operating with real valued spike
    // time may exhibit a different effective refractory time.

    V_.RefractoryCounts_ = nest::Time(nest::Time::ms(P_.TauR_)).get_steps();
    assert(V_.RefractoryCounts_ >= 0);  // since t_ref_ >= 0, this can only fail in error
  }

  /* ----------------------------------------------------------------
   * Update and spike handling functions
   */

  void iaf_psc_alpha2::update(nest::Time const & origin, const nest::long_t from, const nest::long_t to)
  {
    assert(to >= 0 && (nest::delay) from < nest::Scheduler::get_min_delay());
    assert(from < to);

    for ( nest::long_t lag = from ; lag < to ; ++lag )
    {
      if ( S_.r_ == 0 )
      {
        // neuron not refractory
        S_.y3_ = V_.P30_*(S_.y0_ + P_.I_e_)
                 + V_.P31_ex_ * S_.y1_ex_ + V_.P32_ex_ * S_.y2_ex_
                 + V_.P31_in_ * S_.y1_in_ + V_.P32_in_ * S_.y2_in_
                 + V_.expm1_tau_m_ * S_.y3_ + S_.y3_;

        // lower bound of membrane potential
        S_.y3_ = ( S_.y3_ < P_.LowerBound_ ? P_.LowerBound_ : S_.y3_);
      }
      else // neuron is absolute refractory
        --S_.r_;

      // alpha shape EPSCs
      S_.y2_ex_  = V_.P21_ex_ * S_.y1_ex_ + V_.P22_ex_ * S_.y2_ex_;
      S_.y1_ex_ *= V_.P11_ex_;

      // Apply spikes delivered in this step; spikes arriving at T+1 have
      // an immediate effect on the state of the neuron
      V_.weighted_spikes_ex_ = B_.ex_spikes_.get_value(lag);
      S_.y1_ex_ += V_.EPSCInitialValue_ * V_.weighted_spikes_ex_;

      // alpha shape EPSCs
      S_.y2_in_  = V_.P21_in_ * S_.y1_in_ + V_.P22_in_ * S_.y2_in_;
      S_.y1_in_ *= V_.P11_in_;

      // Apply spikes delivered in this step; spikes arriving at T+1 have
      // an immediate effect on the state of the neuron
      V_.weighted_spikes_in_ = B_.in_spikes_.get_value(lag);
      S_.y1_in_ += V_.IPSCInitialValue_ * V_.weighted_spikes_in_;

      // threshold crossing
      if ( S_.y3_ >= P_.Theta_)
      {
        S_.r_  = V_.RefractoryCounts_;
        S_.y3_ = P_.V_reset_;
        // A supra-threshold membrane potential should never be observable.
        // The reset at the time of threshold crossing enables accurate integration
        // independent of the computation step size, see [2,3] for details.

        set_spiketime(nest::Time::step(origin.get_steps()+lag+1));
        nest::SpikeEvent se;
        network()->send(*this, se, lag);
      }

      // set new input current
      S_.y0_ = B_.currents_.get_value(lag);

      // log state data
      B_.logger_.record_data(origin.get_steps() + lag);
    }
  }

  void iaf_psc_alpha2::handle(nest::SpikeEvent& e)
  {
    assert(e.get_delay() > 0);

    const nest::double_t s = e.get_weight() * e.get_multiplicity();

    if(e.get_weight() > 0.0)
      B_.ex_spikes_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), s);
    else
      B_.in_spikes_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), s);
  }

  void iaf_psc_alpha2::handle(nest::CurrentEvent& e)
  {
    assert(e.get_delay() > 0);

    const nest::double_t I = e.get_current();
    const nest::double_t w = e.get_weight();

    B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), w * I);
  }

  void iaf_psc_alpha2::handle(nest::DataLoggingRequest& e)
  {
    B_.logger_.handle(e);
  }

} // namespace
