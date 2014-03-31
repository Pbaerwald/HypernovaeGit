#Calculations the scattering cross section as a function of normalized photon energy
#the energy is normalized to the electron rest mass energy AS SEEN IN THE ELECTRONS REST FRAME
#has some approximations for the non-relativisitic and ultra-relativistic 
#regimes as well as the full KN treatment. 
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
#	E -- Particle kinetic energy  [electron rest mass energy]
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

#funcion which expresses the Maxwell-Juttner distribution 
#the probability distribution for a relativistic gas 
#with Lorentz factors gamma. The imputs are the Lorentz
#factor and the temperature which is scaled for the 
#electron rest mass energy. 

function MJ_distribution(gamma::Float64,T::Float64)
    @assert(gamma >= 1.0);
    @assert(T >= 0.0); 
    beta = sqrt(1.0-1.0/gamma/gamma);
    denom = T*besselk(2,1.0/T);
    denom > 0.0 ? (return gamma*gamma*beta*exp(-gamma/T)/denom) : return 0.0;
end

#returns the single energy emissivity for relativistic Bremsstrahlung
#Inputs:
#       omega -- angular frequency of light [electron rest mass energy] 
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

function get_rel_specific_emissivity(this_omega::Float64,this_theta::Float64)
    max_gamma = maximum([1000.0*this_theta,100.0]); 
    

    denominator, denominator_err= quadgk(a -> MJ_distribution(a,this_theta),1.0,max_gamma);
    if denominator == 0.0
        return 0.0; 
    end
    numerator, numerator_err = quadgk(b -> scattering_cross_section(b)*rel_single_energy_emissivity(this_omega,b)*MJ_distribution(b,this_theta),1.0,max_gamma); 
    
    this_result = numerator/denominator; 
    
    this_result > 0.0 ? (return 9.7411e-42*this_result) : (return 0.0); 
end

#Function that returns the non-relativistic, velocity averaged
#Gaunt factor for a thermal distribution. Different approximations are
#used for different temperature regiemes (See Rybicki Lightman pp 160) 

function get_NR_specific_emissivity(this_nu::Float64,this_theta::Float64,this_temp::Float64)
	#TODO ADD BETTER APPROXIMATION FOR THE GUANT FACTOR
	g_ff = 1.2; 
	return 6.8e-38*g_ff*exp(-this_nu/this_theta)/sqrt(this_temp); 
end


function specific_emissivity(nu::Float64,temp::Float64)
    @assert(nu > 0.0); 
    @assert(temp > 0.0); 
    
    #angular frequency scaled to electron rest mass
    omega = 8.0933e-21*nu; 
    
    #temperature scaled to electron rest mass 
    theta = 1.787387e-10*temp;
    
    theta > 0.01 ? (result = get_rel_specific_emissivity(omega,theta)) : (result = get_NR_specific_emissivity(omega,theta,temp));
    
    @assert(result >= 0.0); 
    
    return result; 
    
end
