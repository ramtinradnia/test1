import math, numpy, random

class SyntheticTypeCurve ():

    def __init__(self, q_zero, q_peak, t_peak, t_plat, k, b, a, p_c):
        # Time (t) - measured in months (m)
        self.t = 0 
        self.t_plat = max(t_plat, 0)

        # Standard Volume Flow (q) - measured in Standard meters cubes per second (Sm^3/s)
        self.q_zero = min(q_zero, q_peak)
        self.q_peak = q_peak
        self.t_peak = max(t_peak, 0)

        # Interpolation degree between linear and sinusoidal ramp to peak (k)
        self.k = min(max(k,0),1)

        # Flow decline exponent (b)
        self.b = min(max(b,0),1)

        # Flow decline constant (a) - a=ln(q_zero / q) / t
        self.a = max(a, 0) 
    
        # Noise parameters
        self.c1 = 0 # std deviation of q(t)
        self.c2 = 0 # sine amplitude, fractional
        self.c3 = 0 # sine frequency 
        self.c4 = 1 # major disturbance average, fractional 
        self.c5 = 0 # major disturbance std deviation, 0 < c5 < c4

        # Noise Probability 
        self.p_dist = 0 # probability of major distubance
        self.p_shut = 0 # probability of shut in 
    
        # Inflow Oerformance Relationship 
        self.p_low = 0 # low pressure, high flow
        self.p_high = 100 # high pressure, low flow
        self.x_low = 1.5 # high flow multiplier 
        self.x_high = 0 # flow flow multiplier 
        
        self.p_ref = 22 #reference pressure, psia
        self.g_rec = 0.85 # gas recoverability 
        self.g_init = 250 # initial gas content 
        self.p_l = 750 # Langmuir pressure, psia
        self.v_l = 550 # Langmuir volume, scf/tonne
        self.n = 0.5 # exponent determining shape of inflow performance relationship
        self.p_c = min(max(k,0),1) # performance coefficient --> make into method to find automatically later


    def build_curve(self, ts):
        # return a list of flows from a list of ts
        qs = []
        for t in ts:
            self.t = t
            qs.append(self.calculate_synthetic_type_curve())
        return qs

    # Section 3

    def calculate_ramp_to_peak_flow(self):
        # Equation 2
        # returns the q_flow at time self.t
        return (self.calculate_linear_ramp_slope() * self.t) + self.q_zero 

    def calculate_linear_ramp_slope(self):
        # Equation 3
        # returns the slope in equation 2 at time self.t
        try:     
            return (self.q_peak - self.q_zero) / self.t_peak
        except Exception as err:
            return 0

    def calculate_sine_ramp_flow(self):
        # Equation 4
        # returns the sine function to peak flow at time self.t
        try:
            return (self.q_peak / 2) * math.sin((math.pi / self.t_peak) * self.t - (math.pi / 2)) + (self.q_peak / 2)
        except Exception as err:
            return 0

    def calculate_ramp_flow(self):
        # Equation 5
        # returns stiffness parameter (k) between two ramping functions at time self.t
        return (self.k * self.calculate_sine_ramp_flow()) + ((1 - self.k) * self.calculate_ramp_to_peak_flow()) 

    def calculate_peak_flow_plateau(self): 
        # Equation 6
        # returns gas flow during plateau flow 
        return (self.q_peak)

    def calculate_decline_flow(self):
        if self.b == 0:
            # Equation 8
            # returns exponential decline flow at time self.t
            return self.q_peak * math.exp((-1 * (self.a)) * (self.t - self.t_peak - self.t_plat))
        else:
            # Equation 7
            # returns decline flow at time self.t 
            return self.q_peak / math.pow((1 + (self.b * self.a * (self.t - self.t_peak - self.t_plat))) , (1/self.b))
                        
    def calculate_synthetic_type_curve(self):
        # Equation 9
        if self.t == 0:
            return self.q_zero
        else:
            if self.t > 0 and self.t <= self.t_peak:
                return self.calculate_ramp_flow()  
            else:
                if self.t > self.t_peak and self.t <= (self.t_peak + self.t_plat):
                    return self.calculate_peak_flow_plateau()    
                else:
                    if self.t > (self.t_peak + self.t_plat):
                        return self.calculate_decline_flow()

    # Section 4

    def calculate_noise(self):
        x_noise = 1

        # Equation 10
        # returns gas flow rate noise
        if self.c1 > 0:
            x_noise = numpy.random.normal(1, self.c1)
        elif self.c2 > 0:
            x_noise = 1.0 + (self.c2 * math.sin(self.c3 * self.t * 2.0 * math.pi))
        return x_noise

    def calculate_disturbance_noise(self):
        random.seed(1000)
        numpy.random.seed(1000)
        x_dist = 1
    
        # Equation 11
        # return gas flow rate disturbance
        if random.random() <= self.p_dist:
            x_dist = numpy.random.normal(self.c4, self.c5)
        return x_dist

    def calculate_shutin_noise(self):
        x_shut = 1
        random.seed(1000)
        numpy.random.seed(1000)

        # Equation 12 
        # returns gas flow rate during shut-in (period of no flow for a set amount of time)
        if random.random() <= self.p_shut:
            x_shut = 0
        return x_shut

    def calculate_synthetic_type_curve_noise(self):
        # Equation 13
        # returns gas flow rate at time self.t with noise, disturbance, and shut-ins 
        return self.calculate_synthetic_type_curve() * self.calculate_noise() * self.calculate_disturbance_noise() * self.calculate_shutin_noise()
                            
    def calculate_IPR(self, q, p): #remove function + variables?
        if p < self.p_low:
            return q * self.x_low
        elif self.p_low <= p <= self.p_ref:
            return q * (1 - ((self.p_ref - p) / (self.p_ref - self.p_low) * (1 - self.x_low)))
        elif self.p_ref <= p <= self.p_high:
            return q * (1 - ((self.p_ref - p) / (self.p_ref - self.p_high) * (1 - self.x_high))) 
        else:
            return q * self.x_high
    
    def get_flow(self, p=None): #remove function?
        if p == None:
            p = self.p_ref
        
        if self.calculate_synthetic_type_curve() < self.q_zero:
            return 0
        else:
            return self.calculate_IPR(self.calculate_synthetic_type_curve_noise(), p)

    # Section 5

    def calculate_gas_flow_cum(self):
        # MMSCF
        q_cum=[]
        for temp in len(range(12 * 30)):
            if temp == 0:
                q_cum[temp] = 0
            else:
                q_cum[temp] = q_cum[temp - 1] + (self.calculate_synthetic_type_curve() * 30 * 0.0000002359705826) # remove manual conversion for engunits later 
        return q_cum[self.t]

        
    def get_gas_flow_cum_max(self):
        # MMSCF
        q_max = 0
        for i in len(range(12 * 30)):
            if q_max < self.calculate_gas_flow_cum():
                q_max = self.calculate_gas_flow_cum()
        return q_max
        
    def calculate_gas_reservoir(self):
        return ((self.get_gas_flow_cum_max() / self.g_rec) - self.calculate_gas_flow_cum())    

    def calculate_gas_content(self):
        # scf/tonnes
        return ((self.calculate_gas_reservoir() * 1000000) / (((self.get_gas_flow_cum_max() * 1000000) / self.g_rec) / self.g_init))

    def calculate_p_reservoir(self):
        # psia
        return (self.p_l / ((self.v_l / self.calculate_gas_content) - 1)) 

    def calculate_gas_flow_rate(self):
        # MMSCFD
        return (self.p_c * (((self.calculate_p_reservoir() ** 2) - (self.p_ref ** 2)) ** self.n)) 

    def get_backpressure(self, interval, flow):
        pass

