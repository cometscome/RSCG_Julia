module rscg
    
    function RSCG(eps,n_omega,left_i,right_j,vec_sigma,A,Ln)

    #--Line 2 in Table III.
        vec_x = zeros(Float64,Ln)
        vec_b =   zeros(Float64,Ln)
        vec_b[right_j] = 1.0
        
        vec_r = zeros(Ln)
        vec_p = zeros(Ln)
        vec_r[right_j] = 1.0
        vec_p[right_j] = 1.0
        alpham = 1.0
        betam = 0.0
    #--
        Sigma = vec_b[left_i] #Line 3. Sigma=V.b. V=v1^T, v1^T=e(j)^T = (0,0,0,0,...,1,...,0,0,0)

        vec_Ap = zeros(Ln)
        vec_g = zeros(Complex{Float64},n_omega)
        vec_rhok = ones(Complex{Float64},n_omega)
        vec_rhokp = ones(Complex{Float64},n_omega)    
        vec_rhokm = ones(Complex{Float64},n_omega)
        vec_alpha = zeros(Complex{Float64},n_omega)
        vec_beta = zeros(Complex{Float64},n_omega)
        vec_Theta = zeros(Complex{Float64},n_omega)
        vec_Pi = ones(Complex{Float64},n_omega)*Sigma
        #---
        #flag = true
        hi = 1.0
    
        ep = 1e-15

        for ite in 1:200
        #while abs(hi) > eps        
            #vec_Ap = 
            A_mul_B!(vec_Ap,A,-vec_p)   
            rsum = vec_r'*vec_r#  dot(vec_r,vec_r) #(rk,rk)
            pAp = vec_p'*vec_Ap
            alpha = rsum/pAp #np.dot(vec_p,vec_Ap) #Line 6 (rk,rk)/(pk,A pk)
#        print alpha,rsum
            
            for i=1:Ln
                vec_x[i] = alpha*vec_p[i] #Line 7
                vec_r[i] += -alpha*vec_Ap[i]#Line 8
            end
#            vec_x += alpha*vec_p #Line 7
#            vec_r += - alpha*vec_Ap #Line 8
            beta = vec_r'*vec_r/rsum #Line9 (r_{k+1},r_{k+1})/(rk,rk)
            vec_p = vec_r + beta*vec_p #Line 10
            Sigma = vec_r[left_i] #Line 11 Sigma=V.r_{k+1}
        
            
        
            #---- Lines 12-17
            for j in 1:n_omega
                update = ifelse(abs(vec_rhok[j]) > ep,true,false)
                if update
                    vec_rhokp[j] = vec_rhok[j]*vec_rhokm[j]*alpham/(vec_rhokm[j]*alpham*(1.0+alpha*vec_sigma[j])+alpha*betam*(vec_rhokm[j]-vec_rhok[j]))#Line 13
                    vec_alpha[j] = alpha*vec_rhokp[j]/vec_rhok[j]#Line 14
                    vec_Theta[j] += vec_alpha[j]*vec_Pi[j] #Line 15
                    vec_beta[j] = ((vec_rhokp[j]/vec_rhok[j])^2)*beta #Line 16
                    vec_Pi[j] = vec_rhokp[j]*Sigma+ vec_beta[j]*vec_Pi[j] #Line 17
                end
                #=
                vec_rhokp[j] = ifelse(update,vec_rhok[j]*vec_rhokm[j]*alpham/(vec_rhokm[j]*alpham*(1.0+alpha*vec_sigma[j])+alpha*betam*(vec_rhokm[j]-vec_rhok[j])), vec_rhok[j])#Line 13
                vec_alpha[j] = ifelse(update, alpha*vec_rhokp[j]/vec_rhok[j],0.0)#Line 14
                vec_Theta[j] = vec_Theta[j]+vec_alpha[j]*vec_Pi[j] #Line 15
                vec_beta[j] = ifelse(update,((vec_rhokp[j]/vec_rhok[j])^2)*beta,1.0) #Line 16
                vec_Pi[j] = ifelse(update,vec_rhokp[j]*Sigma+ vec_beta[j]*vec_Pi[j],vec_Pi[j]) #Line 17
                =#
                vec_g[j] = vec_Theta[j]
                vec_rhokm[j] = vec_rhok[j]
                vec_rhok[j] = vec_rhokp[j]
               
            end
        
            
            #----
            alpham = alpha
            betam = beta
            hi = rsum
            #println(ite," ",vec_r'*vec_r)
            #println(sum(-A*vec_x+vec_b))
 #           println(rsum)


        end
    
        #println(sum(A*vec_x))
        #println("answer? ",sum(A*vec_x+vec_b))
        #println(" d ",vec_x[1])
        return vec_g
    end

    function calc_meanfields_full_finite(A,Nx,Ny,Ln,T,omegamax)

        cc = spzeros(Nx*Ny,Nx*Ny) #lil_matrix((Nx*Ny,Nx*Ny))
#    omegamax = omegac #pi*T(2*n+1), omegac/(T*pi)
        n_omega = Int((Int(omegamax/(T*pi)))/2-1)
    
        vec_sigma = zeros(Complex{Float64},2*n_omega)

        for n=1:2*n_omega
            vec_sigma[n] = π*T*(2.0*(n-n_omega-1)+1)*im
        end
        A = full(A)
        w,v = eig(A)

    
        for ix= 1:Nx#  in range(Nx):
            for iy= 1:Ny#  in range(Ny):
                ii = (iy-1)*Nx+ix                
                jj = ii + Nx*Ny
                right_j = jj
                left_i = ii
                vec_g = calc_green(w,v,2*n_omega,left_i,right_j,vec_sigma,A,Ln)
                vec_g += -1./(vec_sigma.*vec_sigma)
                cc[ii,ii] = real(T*sum(vec_g))-1/(T*4)
                
                #stop
            end
        end

        return cc
    end

    function calc_meanfields_full_finite_HF(A,Nx,Ny,Ln,T,omegamax)

        cdc = spzeros(Nx*Ny,Nx*Ny) #lil_matrix((Nx*Ny,Nx*Ny))
#    omegamax = omegac #pi*T(2*n+1), omegac/(T*pi)
        n_omega = Int((Int(omegamax/(T*pi)))/2-1)
    
        vec_sigma = zeros(Complex{Float64},2*n_omega)

        for n=1:2*n_omega
            vec_sigma[n] = π*T*(2.0*(n-n_omega-1)+1)*im
        end
        A = full(A)
        w,v = eig(A)

    
        for ix= 1:Nx#  in range(Nx):
            for iy= 1:Ny#  in range(Ny):
                ii = (iy-1)*Nx+ix                
                jj = ii# + Nx*Ny
                right_j = jj
                left_i = ii         
                vec_g = calc_green(w,v,2*n_omega,left_i,right_j,vec_sigma,A,Ln)
                vec_g += -1./vec_sigma
                
                cdc[ii,ii] = real(T*sum(vec_g))+1/2
            end
        end

        return cdc
    end

    function calc_green(w,v,n_omega,left_i,right_j,vec_sigma,A,Ln)
        vec_g = zeros(Complex{Float64},n_omega)
        
        for i=1:Ln
            for n=1:n_omega
                vec_g[n] += v[left_i,i]*v[right_j,i]/(vec_sigma[n] - w[i])
            end
        end
        return vec_g
    end

    function calc_meanfields_RSCG_HF(eps,A,Nx,Ny,Ln,T,omegamax)

        cc = spzeros(Nx*Ny,Nx*Ny) #lil_matrix((Nx*Ny,Nx*Ny))
#    omegamax = omegac #pi*T(2*n+1), omegac/(T*pi)
        n_omega = Int((Int(omegamax/(T*pi)))/2-1)
    
        vec_sigma = zeros(Complex{Float64},2*n_omega)
#        shift = 3.0
#        for i=1:Ln
#            A[i,i] = A[i,i] + shift
#        end
    

        for n=1:2*n_omega# in range(2*n_omega):
            #println(π*T*(2.0*(n-n_omega-1)+1)*im)
            vec_sigma[n] = π*T*(2.0*(n-n_omega-1)+1)*im #+shift
        end

    
        for ix= 1:Nx#  in range(Nx):
            for iy= 1:Ny#  in range(Ny):
                ii = (iy-1)*Nx+ix                
                jj = ii #+ Nx*Ny
                right_j = jj#+ Nx*Ny
                left_i = ii#+ Nx*Ny
                vec_g = RSCG(eps,2*n_omega,left_i,right_j,vec_sigma,A,Ln)
                vec_g += -1./vec_sigma
                cc[ii,ii] = real(T*sum(vec_g))+1/2
            end
        end

        return cc
    end
    
            
    function calc_meanfields_RSCG(eps,A,Nx,Ny,Ln,T,omegamax)

        cc = spzeros(Nx*Ny,Nx*Ny) #lil_matrix((Nx*Ny,Nx*Ny))
#    omegamax = omegac #pi*T(2*n+1), omegac/(T*pi)
        n_omega = Int((Int(omegamax/(T*pi)))/2-1)
    
        vec_sigma = zeros(Complex{Float64},2*n_omega)

#        shift = 3.0
#        for i=1:Ln
#            A[i,i] = A[i,i] + shift
#        end
    
        for n=1:2*n_omega# in range(2*n_omega):
            #println(π*T*(2.0*(n-n_omega-1)+1)*im)
            vec_sigma[n] = π*T*(2.0*(n-n_omega-1)+1)*im#+shift
        end
    


    
        for ix= 1:Nx#  in range(Nx):
            for iy= 1:Ny#  in range(Ny):
                ii = (iy-1)*Nx+ix                
                jj = ii + Nx*Ny
                right_j = jj
                left_i = ii
                vec_g = RSCG(eps,2*n_omega,left_i,right_j,vec_sigma,A,Ln)
                vec_g += -1./(vec_sigma.*vec_sigma)
                cc[ii,ii] = real(T*sum(vec_g))-1/(T*4)
            end
        end

        return cc
    end

end