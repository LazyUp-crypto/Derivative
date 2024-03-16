import numpy as np
#from numpy.polynomial import Polynomial

class AmericanOptionsLSMC(object):
    def __init__(self, S0,K, T, M, r, div, sigma, simulations):
        try:
            self.S0 = S0
            assert T>0
            self.T = float(T)
            assert r > 0
            self.r = float(r)
            assert K>0
            self.K = float(K)
            assert M > 0
            self.M = int(M)
            assert sigma > 0
            self.sigma = float(sigma)
            assert div >= 0
            self.div = float(div)
            assert simulations>0
            self.simulations = int(simulations)
        except ValueError:
            print("Error input put options parameters")

        if S0 < 0 or K < 0 or T <= 0 or r < 0 or div < 0 or sigma < 0:
            raise ValueError('Error: Negative inputs are not allowed')
        
        self.time_unit = self.T/float(self.M)
        self.discount = np.exp(- self.r * self.time_unit)

    #标的资产价格stock price path matrix
    @property
    def SPrice_Matrix(self,seed=123):
        np.random.seed(seed)
        #row 0 for S0,row 1 to M+1 for MCPrice,M+1行：时间，simulations 列
        SPrice_Matrix = np.zeros((self.M+1,self.simulations),dtype = np.float64)
        SPrice_Matrix[0,:] = self.S0
        for t in range(1,self.M+1):
            Z = np.random.randn(self.simulations)
            SPrice_Matrix[t,:] = SPrice_Matrix[t-1,:] * np.exp((self.r - 1/2 * self.sigma ** 2) * self.time_unit +  self.sigma * Z * np.sqrt(self.time_unit))
        return SPrice_Matrix
    
    #payoff 结合行权价K,看跌期权
    @property   
    def MCpayoff(self):
        payoff = np.maximum((self.K - self.SPrice_Matrix),np.zeros((self.M+1,self.simulations),dtype=np.float64))
        return payoff
    
    #value
    @property
    def value_vector(self):
        value_matrix = np.zeros_like(self.SPrice_Matrix)
        #从后向前推，最后没有continuation value，所以相当于exercise value
        value_matrix[-1,:] = self.MCpayoff[-1,:]
        for t in range(self.M-1,0,-1):
            #从最后M开始
            #五次多项式拟合，预测在时间t的股票价格和在时间t+1的continuation value（t+1)之间的关系，返回的是多项式的系数
            coeff = np.polynomial.polynomial.polyfit(self.SPrice_Matrix[t,:],value_matrix[t+1,:] * self.discount,5)
            #np.polyval可以计算多项式在给定点的值,coeff是从最低阶到最高阶的系数，所以反转
            #下边计算的就是第t时刻的持续价值
            continuation_value = np.polyval(coeff[::-1],self.SPrice_Matrix[t,:])
            #比较当期行权的payoff和continuation value
            value_matrix[t,:] = np.where(self.MCpayoff[t,:] > continuation_value,
                                         self.MCpayoff[t,:],
                                         value_matrix[t+1,:] * self.discount)
        return value_matrix[1,:] * self.discount
    
    @property
    def price(self):
        return np.sum(self.value_vector) / float(self.simulations)

    @property
    def delta(self):
        diff = self.S0 * 0.01
        myCall_1 = AmericanOptionsLSMC(self.S0 + diff,self.K,self.T,self.M,self.r,self.div,self.sigma,self.simulations)

        myCall_2 = AmericanOptionsLSMC(self.S0 - diff,self.K,self.T,self.M,self.r,self.div,self.sigma,self.simulations)

        return (myCall_1.price - myCall_2.price) / float(2. * diff)
    
    @property
    def gamma(self):
        diff = self.S0 * 0.01
        myCall_1 = AmericanOptionsLSMC(self.S0 + diff, 
                                       self.K, self.T, self.M, 
                                       self.r, self.div, self.sigma, self.simulations)
        myCall_2 = AmericanOptionsLSMC(self.S0 - diff, 
                                       self.K, self.T, self.M, 
                                       self.r, self.div, self.sigma, self.simulations)
        return (myCall_1.delta - myCall_2.delta) / float(2. * diff)
    
    @property 
    def vega(self):
        diff = self.sigma * 0.01
        myCall_1 = AmericanOptionsLSMC(self.S0, 
                                       self.K, self.T, self.M, 
                                       self.r, self.div, self.sigma+diff, self.simulations)
        myCall_2 = AmericanOptionsLSMC(self.S0, 
                                       self.K, self.T, self.M, 
                                       self.r, self.div, self.sigma-diff, self.simulations)
        return (myCall_1.price - myCall_2.price) / float(2. * diff)
    
    @property
    def rho(self):        
        diff = self.r * 0.01
        if (self.r - diff) < 0:        
            myCall_1 = AmericanOptionsLSMC(self.S0, 
                                       self.K, self.T, self.M, 
                                       self.r + diff, self.div, self.sigma, 
                                       self.simulations)
            myCall_2 = AmericanOptionsLSMC(self.S0, 
                                       self.K, self.T, self.M, 
                                       self.r, self.div, self.sigma, 
                                       self.simulations)
            return (myCall_1.price - myCall_2.price) / float(diff)
        else:
            myCall_1 = AmericanOptionsLSMC(self.S0, 
                                       self.K, self.T, self.M, 
                                       self.r + diff, self.div, self.sigma, 
                                       self.simulations)
            myCall_2 = AmericanOptionsLSMC(self.S0, 
                                       self.K, self.T, self.M, 
                                       self.r - diff, self.div, self.sigma, 
                                       self.simulations)
            return (myCall_1.price - myCall_2.price) / float(2. * diff)

    @property
    def theta(self):
        diff = 1/ 252.
        myCall_1 = AmericanOptionsLSMC(self.S0, 
                                       self.K, self.T+diff, self.M, 
                                       self.r, self.div, self.sigma, 
                                       self.simulations)
        myCall_2 = AmericanOptionsLSMC(self.S0, 
                                       self.K, self.T-diff, self.M, 
                                       self.r, self.div, self.sigma, 
                                       self.simulations)
        return (myCall_2.price - myCall_1.price) / float(2. * diff)
    
if __name__ == '__main__':
    AmericanPUT = AmericanOptionsLSMC( 36., 40., 1., 50, 0.06, 0.06, 0.2, 10000)
    print ('Price: ', AmericanPUT.price)
    print ('Delta: ', AmericanPUT.delta)
    print ('Gamma: ', AmericanPUT.gamma)
    print ('Vega:  ', AmericanPUT.vega)
    print ('Rho:   ', AmericanPUT.rho)
    print ('Theta: ', AmericanPUT.theta)