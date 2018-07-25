classdef MPmat
    
    properties (SetAccess = private)
        id
        % matpar : STEEL FIXED PROPERTIES
        Fy;  %  = matpar(1)  : yield stress
        E0;  %  = matpar(2)  : initial stiffness
        b;   %  = matpar(3)  : hardening ratio (Esh/E0)
        R0;  %  = matpar(4)  : exp transition elastic-plastic
        cR1; %  = matpar(5)  : coefficient for changing R0 to R
        cR2; %  = matpar(6)  : coefficient for changing R0 to R
        a1;  %  = matpar(7)  : coefficient for isotropic hardening in compression
        a2;  %  = matpar(8)  : coefficient for isotropic hardening in compression
        a3;  %  = matpar(9)  : coefficient for isotropic hardening in tension
        a4;  %  = matpar(10) : coefficient for isotropic hardening in tension
        sigini;  % initial stress
        % hstvP : STEEL HISTORY VARIABLES
        epsminP; %  = hstvP(1) : max eps in compression
        epsmaxP; %  = hstvP(2) : max eps in tension
        epsplP;  %  = hstvP(3) : plastic excursion
        epss0P;  %  = hstvP(4) : eps at asymptotes intersection
        sigs0P;  %  = hstvP(5) : sig at asymptotes intersection
        epssrP;  %  = hstvP(6) : eps at last inversion point
        sigsrP;  %  = hstvP(7) : sig at last inversion point
        konP;    %  = hstvP(8) : index for loading/unloading
        % hstv : STEEL HISTORY VARIABLES
        epsP;    %  = strain at previous converged step
        sigP;    %  = stress at previous converged step
        eP;      %   stiffness modulus at last converged step;
        
        epsmin;
        epsmax;
        epspl;
        epss0;
        sigs0;
        epsr;
        sigr;
        kon;
        sig;
        e;
        epcs;    %  = strain at current step
    end
    
    methods
        function self = MPmat(Fy_, E0_, b_, R0_, cR1_, cR2_)
            persistent id_
            if (isempty(id_)); id_ = 0; end
            
            self.id = id_;
            id_ = id_ + 1;
            
            self.Fy = Fy_;
            self.E0 = E0_;
            self.b = b_;
            self.R0 = R0_;
            self.cR1 = cR1_;
            self.cR2 = cR2_;
            self.sigini = 0;
            
            self.konP = 0;
            self.kon = 0;
            
            self.a1 = 0;
            self.a2 = 1;
            self.a3 = 0;
            self.a4 = 1;
            
            self.eP = self.E0;
            self.epsP = 0.0;
            self.sigP = 0.0;
            self.sig = 0.0;
            self.epcs = 0.0;
            self.e = self.E0;
            
            self.epsmaxP = self.Fy/self.E0;
            self.epsminP = -self.epsmaxP;
            self.epsplP = 0.0;
            self.epss0P = 0.0;
            self.sigs0P = 0.0;
            self.epssrP = 0.0;
            self.sigsrP = 0.0;
            
        end
        
        function Print(self)
            fprintf('MP material:\n');
            fprintf('fy = %f\n', self.Fy);
            fprintf('E0 = %f\n', self.E0);
            fprintf('b = %f\n', self.b);
            fprintf('R = %f\n', self.R0);
            fprintf('cR1 = %f\n', self.cR1);
            fprintf('cR2 = %f\n', self.cR2);
            fprintf('a1 = %f\n', self.a1);
            fprintf('a2 = %f\n', self.a2);
            fprintf('a3 = %f\n', self.a3);
            fprintf('a4 = %f\n', self.a4);
        end
        
        function Eini = getInitialTangent(self)
            Eini = self.E0;
        end
        
        function sigini = getInitialYS(self)
            sigini = self.Fy;
        end
        
        function bini = getInitialHR(self)
            bini = self.b;
        end
        
        function Rini = getInitialTransEP(self)
            Rini = self.R0;
        end
        
        function self = setParameters(self, E0_, Fy_, b_, R0_)
            self.E0 = E0_;
            self.Fy = Fy_;
            self.b = b_;
            self.R0 = R0_;
        end
        
        
        
        function self = setTrialStrain(self, trialStrain, strainRate)
            Esh = self.b * self.E0;
            epsy = self.Fy / self.E0;
            
            if (self.sigini ~= 0)
                epsini = self.sigini / self.E0;
                self.epcs = trialStrain + epsini;
            else
                self.epcs = trialStrain;
            end
            
            deps = self.epcs - self.epsP;
            
            self.epsmax = self.epsmaxP;
            self.epsmin = self.epsminP;
            self.epspl  = self.epsplP;
            self.epss0  = self.epss0P;
            self.sigs0  = self.sigs0P;
            self.epsr   = self.epssrP;
            self.sigr   = self.sigsrP;
            self.kon    = self.konP;
            
            if (self.kon == 0 || self.kon == 3)
                if (abs(deps) < 10*eps)
                    self.e = self.E0;
                    self.sig = self.sigini;
                    self.kon = 3;
                    return;
                else
                    self.epsmax = epsy;
                    self.epsmin = -epsy;
                    
                    if (deps < 0)
                        self.kon = 2;
                        self.epss0 = self.epsmin;
                        self.sigs0 = -self.Fy;
                        self.epspl = self.epsmin;
                    else
                        self.kon = 1;
                        self.epss0 = self.epsmax;
                        self.sigs0 = self.Fy;
                        self.epspl = self.epsmax;
                    end
                end
            end
            
            if(self.kon == 2 && deps > 0)
                self.kon = 1;
                self.epsr = self.epsP;
                self.sigr = self.sigP;
                % epsmin = min(epsP, epsmin);
                if (self.epsP < self.epsmin)
                    self.epsmin = self.epsP;
                end
                d1 = (self.epsmax - self.epsmin) / (2.0*(self.a4 * epsy));
                shft = 1.0 + self.a3 * d1^0.8;
                self.epss0 = (self.Fy * shft - Esh * epsy * shft - self.sigr + self.E0 * self.epsr) / (self.E0 - Esh);
                self.sigs0 = self.Fy * shft + Esh * (self.epss0 - epsy * shft);
                self.epspl = self.epsmax;
                
            elseif (self.kon == 1 && deps < 0)
                self.kon = 2;
                self.epsr = self.epsP;
                self.sigr = self.sigP;
                %  epsmax = max(epsP, epsmax);
                if (self.epsP > self.epsmax)
                    self.epsmax = self.epsP;
                end
                
                d1 = (self.epsmax - self.epsmin) / (2.0*(self.a2 * epsy));
                shft = 1.0 + self.a1 * d1^0.8;
                self.epss0 = (-self.Fy * shft + Esh * epsy * shft - self.sigr + self.E0 * self.epsr) / (self.E0 - Esh);
                self.sigs0 = -self.Fy * shft + Esh * (self.epss0 + epsy * shft);
                self.epspl = self.epsmin;
                
            end
            
            xi     = abs((self.epspl - self.epss0) / epsy);
            R      = self.R0 * (1.0 - (self.cR1 * xi) / (self.cR2 + xi));
            epsrat = (self.epcs - self.epsr) / (self.epss0 - self.epsr);
            dum1  = 1.0 + abs(epsrat)^R;
            dum2  = dum1^(1/R);
            
            self.sig   = self.b*epsrat + (1.0 - self.b) * epsrat / dum2;
            self.sig   = self.sig * (self.sigs0 - self.sigr) + self.sigr;
            
            self.e = self.b + (1.0 - self.b) / (dum1 * dum2);
            self.e = self.e * (self.sigs0 - self.sigr) / (self.epss0 - self.epsr);
            
        end
        
        function strain = getStrain(self)
            strain = self.epcs;
        end
        
        function stress = getStress(self)
            stress = self.sig;
        end
        
        function Ecs = getTangent(self)
            Ecs = self.e;
        end
        
        function self = commitState(self)
            self.epsminP = self.epsmin;
            self.epsmaxP = self.epsmax;
            self.epsplP = self.epspl;
            self.epss0P = self.epss0;
            self.sigs0P = self.sigs0;
            self.epssrP = self.epsr;
            self.sigsrP = self.sigr;
            self.konP = self.kon;
            
            self.eP = self.e;
            self.sigP = self.sig;
            self.epsP =self.epcs;
        end
        
        function self = perturbE0(self, delta)
            self.E0 = self.E0 + delta;
        end
        
        function self = perturbYS(self, delta)
            self.Fy = self.Fy + delta;
        end
        
        function self = perturbHR(self, delta)
            self.b = self.b + delta;
        end
        
        function self = perturbTransEP(self, delta)
            self.R0 = self.R0 + delta;
        end
        
    end
end