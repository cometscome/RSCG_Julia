
module Rscgsolver
    export calc_meanfields_RSCG_both
    using SparseArrays
    using LinearAlgebra
    using Distributed

    function RSCGs(eps,n_omega,left_i,right_j,vec_sigma,A,Ln)

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
            mul!(vec_Ap, A, -vec_p)
#            A_mul_B!(vec_Ap,A,-vec_p)   
            rsum = vec_r'*vec_r#  dot(vec_r,vec_r) #(rk,rk)
            pAp = vec_p'*vec_Ap
            alpha = rsum/pAp #np.dot(vec_p,vec_Ap) #Line 6 (rk,rk)/(pk,A pk)
#        print alpha,rsum
            
            for i=1:Ln
                vec_x[i] += alpha*vec_p[i] #Line 7
                vec_r[i] += -alpha*vec_Ap[i]#Line 8
            end
#            vec_x += alpha*vec_p #Line 7
#            vec_r += - alpha*vec_Ap #Line 8
            beta = vec_r'*vec_r/rsum #Line9 (r_{k+1},r_{k+1})/(rk,rk)
            vec_p = vec_r + beta*vec_p #Line 10
            Sigma = vec_r[left_i] #Line 11 Sigma=V.r_{k+1}
        
            
        
            #---- Lines 12-17
            for j in 1:n_omega
                update = ifelse(abs(vec_rhok[j]) > eps,true,false)
                if update
                    vec_rhokp[j] = vec_rhok[j]*vec_rhokm[j]*alpham/(vec_rhokm[j]*alpham*(1.0+alpha*vec_sigma[j])+alpha*betam*(vec_rhokm[j]-vec_rhok[j]))#Line 13
                    vec_alpha[j] = alpha*vec_rhokp[j]/vec_rhok[j]#Line 14
                    vec_Theta[j] += vec_alpha[j]*vec_Pi[j] #Line 15
                    vec_beta[j] = ((vec_rhokp[j]/vec_rhok[j])^2)*beta #Line 16
                    vec_Pi[j] = vec_rhokp[j]*Sigma+ vec_beta[j]*vec_Pi[j] #Line 17
                end

                vec_g[j] = vec_Theta[j]
                vec_rhokm[j] = vec_rhok[j]
                vec_rhok[j] = vec_rhokp[j]
               
            end
        
            
            #----
            alpham = alpha
            betam = beta
            hi = rsum


        end
    
        return vec_g
    end

    function RSCG_vec(eps,vec_left,right_j,σ,A,n)
        N=length(σ)
        b = zeros(Float64,n)
        b[right_j] = 1.0

        m = length(vec_left)
        kmax = 20000

    #--Line 2 in Table III.
        x = zeros(Float64,n)
        r = copy(b)  
        p = copy(b)
        αm = 1.0
        βm = 0.0
    #--
        Σ =zeros(Float64,m)
        for mm=1:m
            Σ[mm] = b[vec_left[mm]] #Line 3
        end        
        Θ = zeros(Complex{Float64},N,m)
        Π = ones(Complex{Float64},N,m)
        for mm=1:m
            Π[:,mm] *= Σ[mm]
        end
        ρk = ones(Complex{Float64},N)
        ρkm = copy(ρk)
        ρkp = copy(ρk)
        Ap = similar(p)
    
        for k=0:kmax
            mul!(Ap,A,-p)
            #A_mul_B!(Ap,A,-p)
            rnorm = r'*r
            α = rnorm/(p'*Ap)
            x += α*p #Line 7
            r += -α*Ap #Line 8
            β = r'*r/rnorm #Line9
            p = r + β*p #Line 10
        
            for mm=1:m
                Σ[mm] = r[vec_left[mm]] #Line 11
            end    
            
            for j = 1:N
                update = ifelse(abs(ρk[j]) > eps,true,false)
                if update
                    ρkp[j] = ρk[j]*ρkm[j]*αm/(ρkm[j]*αm*(1.0+α*σ[j])+α*βm*(ρkm[j]-ρk[j]))#Line 13
                    αkj = α*ρkp[j]/ρk[j]#Line 14
                    Θ[j,:] += αkj*Π[j,:] #Line 15
                    βkj = ((ρkp[j]/ρk[j])^2)*β #Line 16
                    Π[j,:] = ρkp[j]*Σ+ βkj*Π[j,:] #Line 17

                end
                ρkm[j] = ρk[j]
                ρk[j] = ρkp[j]
            end
            αm = α
            βm = β
            hi = rnorm
            if hi < eps
                return Θ
            end
        end
    
    
        println("Not converged")
        return Θ

    end



    function calc_meanfields_full_finite(A,Nx,Ny,Ln,T,omegamax)

        cc = spzeros(Nx*Ny,Nx*Ny) #lil_matrix((Nx*Ny,Nx*Ny))
#    omegamax = omegac #pi*T(2*n+1), omegac/(T*pi)
        n_omega = Int((Int(omegamax/(T*pi)))/2-1)
    
        vec_sigma = zeros(Complex{Float64},2*n_omega)

        for n=1:2*n_omega
            vec_sigma[n] = π*T*(2.0*(n-n_omega-1)+1)*im
        end
        A = Matrix(A)
        w,v = eigen(A)

    
        for ix= 1:Nx#  in range(Nx):
            for iy= 1:Ny#  in range(Ny):
                ii = (iy-1)*Nx+ix                
                jj = ii + Nx*Ny
                right_j = jj
                left_i = ii
                vec_g = calc_green(w,v,2*n_omega,left_i,right_j,vec_sigma,A,Ln)
                vec_g += -1 ./(vec_sigma.*vec_sigma)
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
        A = Matrix(A)
        w,v = eigen(A)

    
        for ix= 1:Nx#  in range(Nx):
            for iy= 1:Ny#  in range(Ny):
                ii = (iy-1)*Nx+ix                
                jj = ii# + Nx*Ny
                right_j = jj
                left_i = ii         
                vec_g = calc_green(w,v,2*n_omega,left_i,right_j,vec_sigma,A,Ln)
                vec_g += -1 ./vec_sigma
                
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
                vec_g += -1 ./vec_sigma
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

    
        for n=1:2*n_omega# in range(2*n_omega):
            #println(π*T*(2.0*(n-n_omega-1)+1)*im)
            vec_sigma[n] = π*T*(2.0*(n-n_omega-1)+1)*im#+shift
        end
    
        vec_c = pmap(ii -> calc_meanfields_RSCG_i(ii,vec_sigma,A,Ln,T,Nx,Ny,eps),1:Nx*Ny)
        for ii=1:Nx*Ny
            cc[ii,ii] = vec_c[ii]
            #calc_meanfields_RSCG_i(ii,vec_sigma,A,Ln,T,Nx,Ny,eps)
        end

        #=
        for ix= 1:Nx#  in range(Nx):
            for iy= 1:Ny#  in range(Ny):
                ii = (iy-1)*Nx+ix                
                jj = ii + Nx*Ny
                right_j = jj
                left_i = ii
                vec_left = [left_i]
                vec_g = RSCG_vec(eps,vec_left,right_j,vec_sigma,A,Ln)[:,1]
            
#                vec_g = RSCG(eps,2*n_omega,left_i,right_j,vec_sigma,A,Ln)
                vec_g += -1 ./(vec_sigma.*vec_sigma)
                cc[ii,ii] = real(T*sum(vec_g))-1/(T*4)
            end
        end
    =#

        return cc
    end

    function calc_meanfields_RSCG_i(ii,vec_sigma,A,Ln,T,Nx,Ny,eps)
        jj = ii + Nx*Ny
        right_j = jj
        left_i = ii
        vec_left = [left_i]
        vec_g = RSCG_vec(eps,vec_left,right_j,vec_sigma,A,Ln)[:,1]
            
#                vec_g = RSCG(eps,2*n_omega,left_i,right_j,vec_sigma,A,Ln)
        vec_g += -1 ./(vec_sigma.*vec_sigma)
        return real(T*sum(vec_g))-1/(T*4)
    end

    

    function calc_meanfields_RSCG_both(eps,A,Nx,Ny,Ln,T,omegamax)

        cc = spzeros(Nx*Ny,Nx*Ny) #lil_matrix((Nx*Ny,Nx*Ny))
        cdc = spzeros(Nx*Ny,Nx*Ny)
#    omegamax = omegac #pi*T(2*n+1), omegac/(T*pi)
        n_omega = Int((Int(omegamax/(T*pi)))/2-1)
    
        vec_sigma = zeros(Complex{Float64},2*n_omega)

        for n=1:2*n_omega# in range(2*n_omega):
            #println(π*T*(2.0*(n-n_omega-1)+1)*im)
            vec_sigma[n] = π*T*(2.0*(n-n_omega-1)+1)*im#+shift
        end
    
        vec_cs = pmap(ii -> calc_meanfields_RSCG_both_i(ii,vec_sigma,A,Ln,T,Nx,Ny,eps),1:Nx*Ny)
    
        for ii=1:Nx*Ny
            cc_i,cdc_i = vec_cs[ii]#calc_meanfields_RSCG_both_i(ii,vec_sigma,A,Ln,T,Nx,Ny,eps)        
            cc[ii,ii] = cc_i
            cdc[ii,ii] = cdc_i
        end

        #=
        for ix= 1:Nx#  in range(Nx):
            for iy= 1:Ny#  in range(Ny):
                ii = (iy-1)*Nx+ix                
                jj = ii + Nx*Ny
                right_j = jj
                left_i = ii
                vec_left = [ii,ii+Nx*Ny]
                vec_g = RSCG_vec(eps,vec_left,right_j,vec_sigma,A,Ln)

                vec_g[:,1] += -1 ./(vec_sigma.*vec_sigma)
                vec_g[:,2] += -1 ./(vec_sigma)
                cc[ii,ii] = real(T*sum(vec_g[:,1]))-1/(T*4)
                cdc[ii,ii] = real(T*sum(vec_g[:,2]))+1/2
                cdc[ii,ii] = 1-cdc[ii,ii] 
            end
        end
        =#
        return cc,cdc
    end



    function calc_meanfields_RSCG_both_i(ii,vec_sigma,A,Ln,T,Nx,Ny,eps)
        jj = ii + Nx*Ny
        right_j = jj
        left_i = ii
        vec_left = [ii,ii+Nx*Ny]
        vec_g = RSCG_vec(eps,vec_left,right_j,vec_sigma,A,Ln)

        vec_g[:,1] += -1 ./(vec_sigma.*vec_sigma)
        vec_g[:,2] += -1 ./(vec_sigma)
        cc_i = real(T*sum(vec_g[:,1]))-1/(T*4)
        cdc_i = real(T*sum(vec_g[:,2]))+1/2
        cdc_i = 1-cdc_i
        return cc_i,cdc_i
    end


end