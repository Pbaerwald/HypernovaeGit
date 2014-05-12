#include some physical constants 
const c = 2.99792458e10		#speed of light [cm/s]
const h = 6.6260755e-27		#Planck's constant [erg s]
const k_B = 1.380658e-16	#Boltzmann's constant [erg K^-1] 
const m_e = 9.1093897e-28 	#electron mass [g] 
const charge_e = 4.8032068e-10 	#electron charge [g^1/2 cm^3/2 s^-1 OR erg^1/2 cm^1/2] 

#Calculates the scattering cross section as a function of photon energy
#the energy is normalized to the electron rest mass energy AS SEEN IN THE ELECTRONS REST FRAME
#has some approximations for the non-relativisitic and ultra-relativistic 
#regimes as well as the full KN treatment. 
#
#Inputs: 
#   sx -- photon energy seen in the electrons rest frame in units of electron rest mass energy [1] 
#
#Outputs: 
#   Unnormalized photon-electron scattering cross section [1]
function scattering_cross_section(sx::Float64)
    @assert(sx > 0.0); 
    low_limit = 0.1; 
    high_limit = 100.0; 
    result = 0.0; 
    
    if(sx < low_limit)
       result = 1.0-2*sx+26.0*sx*sx/5; 
    elseif(sx > high_limit)
        result = 3*(log(2*sx)+0.5)/(8*sx); 
    else
        temp0 = log(1.0+2.0*sx);
        temp1 = 2*sx*(1+sx)/(1+2*sx) - temp0; 
        temp1 = (1.0+sx)*temp1/sx/sx/sx; 
        
        temp2 = 0.5*temp0/sx;
            
        temp3 = (1.0+3.0*sx)/(1.0+2.0*sx)/(1.0+2.0*sx); 
            
        result = 0.75(temp1+temp2-temp3);
    end
        
    #TODO add the normalization of the Thompson cross section 
    return result;
end

#function which describes the emission by a single particle 
#function of impact parameter (b), lorentz factor (gamma), angular frequency of emitted light (omega)
function rel_single_particle_emission(x::Float64)
    @assert(x > 0.0)
    bessel = besselk(1,x); 
    return x*bessel*bessel; 
end

#function which expresses the Maxwell-Botlzmann distribution 
#the particle energy probability distribution for a non-relativistic gas 
#as a function of the gas temperatue
#Inputs: 
#	E -- Particle kinetic energy  [m_e c^2]
#	T -- gas temperatue  [electron rest mass energy]
#
#Outputs: 
#	unnormalized Boltzman probability of a particle with energy E in a gas
#	of temperature T 
function MB_distribution(E::Float64,T::Float64)
	@assert(E >= 0); 
	@assert(T >= 0); 
	
	result = 2.0*sqrt(E/pi)*T^(-3/2)*exp(-E/T); 
	
	result >= 0.0 ? (return result) : (return 0.0);
end

function MB_vel_distribution(this_v::Float64, this_T::Float64)
	@assert(this_v > 0.0); 
	@assert(this_T > 0.0); 
	
	this_x = -0.5*m_e*this_v*this_v/(k_B*this_T); 
	
	return this_v*this_v*exp(this_x); 
end
	
	

#funcion which expresses the Maxwell-Juttner distribution 
#the probability distribution for a relativistic gas 
#
#Inputs:
#   gamma -- particles lorentz factor [1] 
#   T     -- temperature of relativistic gas [m_e c^2]
#
#Outputs: 
#   Unnormalized probability of a finding a particle of energy with lorentz factor gamma in a gas of temperature T
function MJ_distribution(gamma::Float64,T::Float64)
    @assert(gamma >= 1.0);
    @assert(T >= 0.0); 
    beta = sqrt(1.0-1.0/gamma/gamma);
    denom = T*besselk(2,1.0/T);
    denom > 0.0 ? (return gamma*gamma*beta*exp(-gamma/T)/denom) : return 0.0;
end

#returns the single energy emissivity for relativistic Bremsstrahlung
#Inputs:
#      omega -- angular frequency of light [m_e c^2]
#
#      gamma -- lorentz factor of electron [1] 
#Ouputs: 
#       gives dW/dt_dV_domega (omega,gamma) from notes 
# NOT YET NORMALIZED 
function rel_single_energy_emissivity(this_omega::Float64,this_gamma::Float64,tol::Float64 = 1e-5)
    @assert(this_omega > 0.0); 
    @assert(this_gamma > 1.0); 
    
    x_min = this_omega/this_gamma/this_gamma;
    
    #needs to be "close" to infinity 
    xf = 10;     
    
    tol = 1e-5; 
    result1, err1= quadgk(rel_single_particle_emission,x_min,0.75*xf);
    result2, err2 = quadgk(rel_single_particle_emission,0.75*xf,xf);
    #make sure the upper bound is close to infinity 
    @assert(result2 < tol) 
    @assert(maximum([err1,err2]) < tol);
    
    return result1+result2;  
end


#returns the bremsstralung specific emissivity given a relativistic thermal plasma
#
#Inputs: 
#   this_omega -- angular frequency of light emitted [m_e c^2] 
#   this_theta -- temperature of plasma [m_e c^2] 
#
#Outputs: 
#   dW/dt_dV_domega integrated over a thermal distribution (omega,T) see notes
#
function get_rel_specific_emissivity(this_omega::Float64,this_theta::Float64)
    max_gamma = maximum([1000.0*this_theta,100.0]); 
    #rel_normalization = 6.1206e-43; 	#normalization for ergs s^-1 cm^-3 Hz^-1 
    rel_normalization = 9.2371e-17; 	#normalization for GeV s^-1 cm^-3 GeV^-1
    

    denominator, denominator_err= quadgk(a -> MJ_distribution(a,this_theta),1.0,max_gamma);
    if denominator == 0.0
        return 0.0; 
    end
    numerator, numerator_err = quadgk(b -> scattering_cross_section(this_omega/b)*rel_single_energy_emissivity(this_omega,b)*MJ_distribution(b,this_theta),1.0,max_gamma); 
    
    this_result = numerator/denominator; 
    
    this_result > 0.0 ? (return rel_normalization*this_result) : (return 0.0); 
end

#returns the more exact bremsstralung specific emissivity geven a non-relativistic thermal plasma
#
#Inputs:
#
#   this_nu     -- frequency of light emitted [m_e c^2]
#   this_theta  -- temperature of plamsa [m_e c^2]
#   this_temp   -- temperature of plasma [K]
#
#Outputs 
#   dW/dt_dV_dnu integrated over a thermal distribution (nu,T) 
#
#Normalization factor:
#FUNCTION YOU ARE WORKING ON 
function get_NR_specific_emissivity(this_nu::Float64,this_temp::Float64)
	NR_normalization = 1.8362e-32;	#normalization for ergs s^-1 cm^-3 Hz^-1
	NR_normalization = 2.7712e-6; 	#normalization for GeV s^-1 cm^-3 GeV^-1 
	
	min_v = sqrt(2.0*this_nu)*3e8; 
	max_v = (100.0*min_v > 0.1*c) ? 0.1*c : 100.*min_v; 
	crossover_v = 4*charge_e^2/(pi*h); 
	
	denominator, denominator_err = quadgk(a -> MB_distribution(a,this_temp),0.0,max_v); 
	if denominator <= 0.0 
		return 0.0;
	end
	
	numerator, numerator_err = quadgk(b -> 
end



#returns the approximate bremsstralung specific emissivity geven a non-relativistic thermal plasma
#
#Inputs:
#
#   this_nu     -- frequency of light emitted [m_e c^2]
#   this_theta  -- temperature of plamsa [m_e c^2]
#   this_temp   -- temperature of plasma [K]
#
#Outputs 
#   dW/dt_dV_dnu integrated over a thermal distribution (nu,T) 
#
#Normalization factor:

function get_approx_NR_specific_emissivity(this_nu::Float64,this_theta::Float64,this_temp::Float64)
	#TODO ADD BETTER APPROXIMATION FOR THE GUANT FACTOR
	g_ff = 1.2; 
	return 6.8e-38*g_ff*exp(-this_nu/this_theta)/sqrt(this_temp); 
end

#returns normalized bremsstralung specific emissivity of a plasma
#determines which approximation to use (NR or relativistic) 
#
#Inputs:
#   nu      -- frequency of light emitted [Hz]
#   temp    -- temperature of plasma considered [K] 
#
#Output: 
#   normalized specific emissivity [erg s^-1 cm^-3 Hz^-1]
#   note that the pre-factor Z^2 n_e n_i is NOT included
#
#   Z       -- ion charge [1]
#   n_e     -- electron density in plasma frame [cm^-3]
#   n_i     -- ion density in plasma frame [cm^-3]


function specific_emissivity(nu::Float64,temp::Float64)
    @assert(nu > 0.0); 
    @assert(temp > 0.0); 
    
    #angular frequency scaled to electron rest mass
    #numerical factor h/(m_e c^2):
    #h      -- Planck's constant    [erg s]
    #m_e    -- electron rest mass   [g]
    #c      -- speed of light       [cm s^-1] 

    omega = 8.0933e-21*nu; 
    
    #temperature scaled to electron rest mass 
    #numerical factor k_b/(m_e c^2):
    #k_b    -- Boltzmann's constant [erg K^-1]
    theta = 1.787387e-10*temp;
    
    theta > 0.01 ? (result = get_rel_specific_emissivity(omega,theta)) : (result = get_approx_NR_specific_emissivity(omega,temp));
    
    @assert(result >= 0.0); 
    
    return result; 
    
end
