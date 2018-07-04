include("ChebyshevPolynomial.jl")
include("RSCG.jl")


module ChevplusRSCG
    import chebyshev
    import rscg
    import Plots

    function calc_A(Nx,Ny,μ,Δ,aa)
        Ln = Nx*Ny*2
        A = spzeros(Ln,Ln)
    
        for ix=1:Nx
            for iy=1:Ny
                #Diagonal element
                ii = (iy-1)*Nx+ix
                jx = ix
                jy = iy
                jj = (jy-1)*Nx+jx
                A[ii,jj] = -μ
                #+1 in x direction
                jx = ifelse(ix == Nx,1,ix+1)
                jy = iy
                jj = (jy-1)*Nx+jx
                A[ii,jj] = -1.0
                #-1 in x direction
                jx = ifelse(ix == 1,Nx,ix-1)
                jy = iy
                jj = (jy-1)*Nx+jx
                A[ii,jj] = -1.0
                #+1 in y direction
                jx = ix
                jy = ifelse(iy == Ny,1,iy+1)
                jj = (jy-1)*Nx+jx
                A[ii,jj] = -1.0
                #-1 in y direction
                jx = ix
                jy = ifelse(iy == 1,Ny,iy-1)
                jj = (jy-1)*Nx+jx
                A[ii,jj] = -1.0            
                            
            end
        end
    
        for ii=1:Nx*Ny
            for jj=1:Nx*Ny
                A[ii+Nx*Ny,jj+Nx*Ny] = -conj(A[ii,jj])
                A[ii,jj+Nx*Ny] = Δ[ii,jj]
                A[ii+Nx*Ny,jj] = conj(Δ[jj,ii])
            end
        end
    
        return A/aa
    
    end

    function update_A!(A,Nx,Ny,μ,Δ,aa)
        for ii=1:Nx*Ny
            for jj=1:Nx*Ny
                A[ii,jj+Nx*Ny] = Δ[ii,jj]/aa
                A[ii+Nx*Ny,jj] = conj(Δ[jj,ii])/aa
            end
        end
    end

    
    function iteration(nc,Nx,Ny,aa,bb,ωc,U,initialΔ,μ,full,RSCG,finite,T,omegamax,itemax)

        Δ = speye(Nx*Ny,Nx*Ny)*initialΔ
        Δold = speye(Nx*Ny,Nx*Ny)*initialΔ     
        A = calc_A(Nx,Ny,μ,Δ,aa)        

        mat_Δ = zeros(typeof(Δ[1,1]),Nx,Ny)    
        
        for ite=1:itemax

        
            if full
                if finite
                    #println("Green's function based full diazonalization")
                    @time Δ = rscg.calc_meanfields_full_finite(aa*A,Nx,Ny,Nx*Ny*2,T,omegamax) 
                else
                    #println("Full diagonalization")
                    @time Δ = chebyshev.calc_meanfields(A,Nx,Ny,ωc) 
                end
            else
                if RSCG
                    #println("RSCG")
                    @time Δ =rscg.calc_meanfields_RSCG(1e-15,aa*A,Nx,Ny,Nx*Ny*2,T,omegamax)
                else
                    #println("Chebyshev")
                   @time Δ = chebyshev.calc_meanfields(nc,A,Nx,Ny,aa,bb,ωc)
                end
            end
                      
            Δ = Δ*U
            update_A!(A,Nx,Ny,μ,Δ,aa)
        
            eps = 0.0
            nor = 0.0
            for i=1:Nx*Ny
                eps += abs(Δ[i,i]-Δold[i,i])^2
                nor += abs(Δold[i,i])^2
            end
            eps = eps/nor
            println("ite = ",ite," eps = ",eps)
            if eps <= 1e-6
                println("End ",Δ[div(Nx,2),div(Ny,2)])
                break
            end
            Δold = Δ
        
            fp = open("./gap.dat","w")
            for ix=1:Nx
                for iy=1:Ny
                    ii = (iy-1)*Nx + ix
                    mat_Δ[ix,iy] = Δ[ii,ii]
                    println(fp,ix,"\t",iy,"\t",abs(mat_Δ[ix,iy]),"\t",angle(mat_Δ[ix,iy]))  
                end
                 println(fp," ")
            end
            println("center (Nx/2,Ny/2): ",mat_Δ[div(Nx,2),div(Ny,2)])
            println("corner (1,1): ",mat_Δ[1,1])
            #plot(mat_Δ)
        
             
            close(fp)
        
        
        end


    
        return mat_Δ
    
    end



end