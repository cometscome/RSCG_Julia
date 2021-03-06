
module Chebyshev
    export calc_meanfields,calc_meanfields_HF
    using SparseArrays
    using LinearAlgebra
    using Distributed



    function calc_meanfields(nc,A,Nx,Ny,aa,bb,ωc)
        cc = spzeros(typeof(A[1,Nx*Ny]),Nx*Ny,Nx*Ny)
        dj = Nx*Ny
        vec_c = pmap(ii -> calc_meanfield_i(ii,dj,nc,A,aa,bb,ωc,Nx,Ny),1:Nx*Ny)
        for ii=1:Nx*Ny
            cc[ii,ii] = vec_c[ii]
        #calc_meanfield_i(ii,nc,A,aa,bb,ωc,Nx,Ny)
        end
        #=
        for ix=1:Nx
            for iy=1:Ny
                ii = (iy-1)*Nx+ix
                jj = ii + Nx*Ny
                vec_ai = calc_polynomials(nc,ii,jj,A)
                cc[ii,ii] = calc_meanfield(vec_ai,aa,bb,ωc,nc)
            end
        end
    =#
        return cc
    end

    function calc_meanfield_i(ii,dj,nc,A,aa,bb,ωc,Nx,Ny)
        jj = ii + dj
        vec_ai = calc_polynomials(nc,ii,jj,A)
        return calc_meanfield(vec_ai,aa,bb,ωc,nc)
    end

    function calc_meanfields_HF(nc,A,Nx,Ny,aa,bb,ωc)
        cdc = spzeros(typeof(A[1,Nx*Ny]),Nx*Ny,Nx*Ny)
        dj = 0
        vec_c = pmap(ii -> calc_meanfield_i(ii,dj,nc,A,aa,bb,ωc,Nx,Ny),1:Nx*Ny)
        for ii=1:Nx*Ny
            cdc[ii,ii] = vec_c[ii]
        #calc_meanfield_i(ii,nc,A,aa,bb,ωc,Nx,Ny)
        end
        #=
        for ix=1:Nx
            for iy=1:Ny
                ii = (iy-1)*Nx+ix
                jj = ii# + Nx*Ny
                vec_ai = calc_polynomials(nc,ii,jj,A)
                cdc[ii,ii] = calc_meanfield(vec_ai,aa,bb,ωc,nc)
            end
        end
    =#
        return cdc
    end

    function calc_meanfields(A,Nx,Ny,ωc)
       
        cc = spzeros(typeof(A[1,Nx*Ny]),Nx*Ny,Nx*Ny)
        A = Matrix(A)
        w,v = eigen(A)
        
        
    
        for ix=1:Nx
            for iy=1:Ny
                ii = (iy-1)*Nx+ix
                jj = ii + Nx*Ny
            
                cc[ii,ii] = 0.0
                for i=1:Nx*Ny*2
                    if w[i] <= 0.0
                        if abs(w[i]) <= ωc
                            cc[ii,ii]+=conj(v[ii,i])*v[jj,i]
                        end
                    end
                end
            end
        end
        A = sparse(A)
        return cc
    end

    function calc_meanfields_HF(A,Nx,Ny,ωc)
       
        cdc = spzeros(typeof(A[1,Nx*Ny]),Nx*Ny,Nx*Ny)
        A = Matrix(A)
        w,v = eigen(A)
        
        
    
        for ix=1:Nx
            for iy=1:Ny
                ii = (iy-1)*Nx+ix
                jj = ii# + Nx*Ny
            
                cdc[ii,ii] = 0.0
                for i=1:Nx*Ny*2
                    if w[i] <= 0.0
                        if abs(w[i]) <= ωc
                            cdc[ii,ii]+=conj(v[ii,i])*v[jj,i]
                        end
                    end
                end
            end
        end
        A = sparse(A)
        return cdc
    end

    function calc_meanfield(vec_ai,aa,bb,ωc,nc)
        ba = acos(-bb/aa)
        ωb = acos(-(ωc+bb)/aa)
        density = vec_ai[0+1]*(ωb-ba)/2
        for nn=1:nc-1
            density += vec_ai[nn+1]*(sin(nn*ωb)-sin(nn*ba))/nn
        end
        density = density*2/π
    end

    function calc_polynomials(nc,left_i,right_j,A)
        Ln = size(A,1) 
        typeA = typeof(A[1,1])
        vec_jnmm = zeros(typeA,Ln)
        vec_jnm = zeros(typeA,Ln)
        vec_jn = zeros(typeA,Ln)
        vec_jn[right_j] = 1.0
        vec_ai = zeros(typeA,nc)
        @inbounds for nn=0:nc-1
            if nn == 0
                vec_jn[right_j] = 1.0
            elseif nn == 1
                mul!(vec_jn, A, vec_jnm)
                #A_mul_B!(vec_jn,A,vec_jnm)
            else
                mul!(vec_jn, A, vec_jnm)
                #A_mul_B!(vec_jn,A,vec_jnm)
#                vec_jn *= 2
                for i=1:Ln
                    vec_jn[i] = vec_jn[i]*2 -vec_jnmm[i]
                end                                               
            end
             
            vec_ai[nn+1] = vec_jn[left_i]
            #println("vec_a, ",vec_ai[nn+1])
            for i=1:Ln
                vec_jnmm[i] = vec_jnm[i]
                vec_jnm[i] = vec_jn[i]
            end


        end
        return vec_ai
    end
end