import matplotlib.pyplot as plt
import numpy as np

class point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    
    def __add__(self, other):
        return point(self.x + other.x, self.y + other.y)

    def __sub__(self, other):
        return point(self.x - other.x, self.y - other.y)

    def print(self):
        print(self.x, self.y)

class instance:
    def __init__(self, n, leng, ang, point):
        self.n = n
        self.leng = leng
        self.ang = ang
        self.pert = np.zeros(n)
        self.point = point
        self.max = sum(leng)

    def draw(self):
        curr = point(0,0)
        ang_sum = 0.0
        for i in range(self.n):
            ang_sum += self.ang[i]
            next = point(curr.x + self.leng[i]*np.cos(ang_sum),curr.y + self.leng[i]*np.sin(ang_sum))
            plt.plot([curr.x, next.x], [curr.y, next.y])
            curr = next
        plt.plot(self.point.x, self.point.y, 'ro')
        plt.xlim(-self.max, self.max)
        plt.ylim(-self.max, self.max)

        plt.show()

    def F(self, pert = False):
        ang_sum = 0.0
        if pert:
            ang = self.ang + self.pert
        else:
            ang = self.ang
        result = point(0,0)
        for i in range(self.n):
            ang_sum += ang[i]
            result.x += self.leng[i] * np.cos(ang_sum)
            result.y += self.leng[i] * np.sin(ang_sum)
        return result

    def G(self, pert = False):
        F_val = self.F(pert = pert)
        return 0.5 * ((F_val.x - self.point.x)**2 + (F_val.y - self.point.y)**2)  
    def getDerivatives(self, pert = False):
        gradG_x = np.zeros(self.n)
        gradG_y = np.zeros(self.n)
        gradG = np.zeros(self.n)
        hessG_x = np.zeros((self.n,self.n))
        hessG_y = np.zeros((self.n,self.n))
        hessG = np.zeros((self.n,self.n))
        if pert:
            ang = self.ang + self.pert
        else:
            ang = self.ang
        ang_sum = 0
        F_x = 0; F_y = 0
        for i in range(self.n):
            ang_sum += ang[i]
            temp_x = self.leng[i]*np.cos(ang_sum)
            temp_y = self.leng[i]*np.sin(ang_sum)
            F_x += temp_x
            F_y += temp_y
            for j in range(i+1):
                gradG_x[j] -= temp_y
                gradG_y[j] += temp_x
                hessG_x[i,j] -= temp_x
                hessG_x[j,i] -= temp_x
                hessG_y[i,j] -= temp_y
                hessG_y[j,i] -= temp_y
        gradG = (F_x - self.point.x)* gradG_x + (F_y - self.point.y)*gradG_y
        hessG = (F_x - self.point.x)* hessG_x + np.outer(gradG_x,gradG_x) + (F_y - self.point.y)*hessG_y + np.outer(gradG_y,gradG_y)
        return 0.5 * ((F_x - self.point.x)**2 + (F_y - self.point.y)**2), gradG, hessG
    
    def steepest(self):
        MAXITER= 100
        TOL = 1e-10
        c = 0.9
        rho = 0.5
        for i in range(MAXITER):
            G, gradG, hessG = self.getDerivatives()
            if (np.linalg.norm(gradG) < TOL):
                break
            alpha = 1.0
            p = - alpha * gradG
            self.pert = p
            G_n = self.G(pert = True)
            k = 0
            while (G_n - G > - c*alpha*np.linalg.norm(p)**2 and k < MAXITER):
                alpha = rho*alpha
                p = - alpha * gradG
                self.pert = p
                G_n = self.G(pert = True)
                k+=1
            self.ang = self.ang + alpha*p
    
    def FR(self):
        MAXITER= 1000
        TOL = 1e-10
        c = 0.9
        rho = 0.8
        G, gradG, hessG = self.getDerivatives()
        p = -gradG
        for i in range(MAXITER):
            if (np.linalg.norm(gradG) < TOL):
                break
            alpha = 1.0
            self.pert = alpha * p
            G_n = self.G(pert = True)
            k = 0
            while (G_n - G > - c*alpha*np.linalg.norm(p)**2 and k < MAXITER):
                alpha_prev = alpha
                alpha = rho*alpha
                p = - alpha/alpha_prev * p  
                self.pert =  p
                G_n = self.G(pert = True)
                k +=1
            print(k)
            self.ang = self.ang + alpha*p
            G_n, gradG_n, _ = self.getDerivatives()
            beta = np.inner(gradG_n,gradG_n)/np.inner(gradG,gradG)
            p = - gradG_n + beta*p
            G = G_n
            gradG = gradG_n
    def newton(self):
        MAXITER = 100
        TOL = 1e-10
        c = 0.9
        rho = 0.8
        G, gradG, hessG = self.getDerivatives()
        p = -gradG
        for i in range(MAXITER):
            if (np.linalg.norm(gradG) < TOL):
                break
            alpha = 1.0
            self.pert = alpha * p
            G_n = self.G(pert = True)
            k = 0
            while (G_n - G > - c*alpha*np.linalg.norm(p)**2 and k < MAXITER):
                alpha_prev = alpha
                alpha = rho*alpha
                p = - alpha/alpha_prev * p  
                self.pert =  p
                G_n = self.G(pert = True)
                k +=1
            print(k, alpha)
            self.ang = self.ang + alpha*p
            G, gradG, hessG = self.getDerivatives()
            p = np.linalg.solve(hessG, - gradG)
            #print(G)
            print(gradG)
ex = instance(3,[3,2,2], [0,0,0], point(0,0.0))
ex.draw()
#ex.G().print()
ex.steepest()
ex.draw()
