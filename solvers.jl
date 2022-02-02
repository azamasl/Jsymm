function lin_search(x,y,obj, F, normF, dir, c1)
    n = length(x)
    m = length(y)
    beta = 0.5
    inner_max_it = 30
    t = 1
    dir_x = dir[1:n]
    dir_y = dir[n+1:n+m]
    x_p = x +t*dir_x
    y_p = y +t*dir_y
    F_p = [obj.∇xL(x_p,y_p); -obj.∇yL(x_p,y_p)]
    normF_p = LinearAlgebra.norm(F_p)
    RHS = c1*normF
    #println("Right hand side is =$RHS")
    in_it = 0

    while  normF -  normF_p < RHS && in_it < inner_max_it
        t = beta*t
        x_p = x +t*dir_x
        y_p = y +t*dir_y
        F_p = [obj.∇xL(x_p,y_p); -obj.∇yL(x_p,y_p)]
        normF_p = LinearAlgebra.norm(F_p)
        in_it = in_it + 1
    end
    return t
end

# Computes the second direction length in the dogleg direction
function getAlpha(p,q,Del)
    #"p'p + alpha^2q'q + 2alpha p'q = Del^2  "
    pdotq = sum(p.*q)
    aa = sum(q.*q)
    bb = 2*pdotq
    cc = sum(p.*p)-Del^2
    alpha = (-bb + sqrt(bb^2 -4*aa*cc))/(2*aa)
    return alpha
end

#######################################
# EGM with fixed stepsize: assuming that nablaF is constant
########################################
function EGM(x,y,obj,dummy0,sp,max_it,prt, use_ls,Ftol)
    # nablaF = [sp.B  sp.A'
    #         -sp.A    sp.C]
    println("dummy print, to make sure next line is not being eatten up")
    println(" Extra Gradient Method with optimal stepsize")
    nablaF  = obj.∇F(x,y)
    eta  = 1/LinearAlgebra.norm(nablaF)
    val, ngx, ngy =0,0,0
    iter=0
    F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
    normF=LinearAlgebra.norm(F)
    normFAll=[]
    Zs=[x;y]
    while ( normF> Ftol ) && iter < max_it
        xk = x -eta*obj.∇xL(x,y)
        yk = y +eta*obj.∇yL(x,y)

        x= x -eta*obj.∇xL(xk,yk)
        y= y +eta*obj.∇yL(xk,yk)

        #nabla_F  = obj.∇F(x,y)
        #eta  = 1/LinearAlgebra.norm(nabla_F)

        Zs = hcat(Zs,[x;y])
        F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
        append!(normFAll, normF)
        normF = LinearAlgebra.norm(F)

        val =obj.L(x,y)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y))
        iter=iter+1
        if prt==1
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            println("iter $iter")
        end
    end
    ng = sqrt(ngx^2+ngy^2)
    return x,y,iter,normFAll, val, ng,Zs
end

#######################################
# Good Broyden method
#######################################
function Broyden(x,y,obj,fixed_stepsz,sp,itNum, prt,use_ls, Ftol)
    println("Broyden's method with using line-search = $use_ls")
    val, ngx, ngy =0,0,0
    n = length(x);m = length(y);
    H = Diagonal(rand(n+m))
    #H = I(n+m)
    F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
    F_old = F
    normF=LinearAlgebra.norm(F)
    normFAll=[];z = [x;y];Zs=z
    iter=0
    while ( normF> Ftol ) && iter < itNum
        z = [x;y]
        p = -H*F
        if(use_ls ==1) # line-search is being used
            gam = lin_search(x,y,obj,F,normF,p,c1)
            if gam >= stepsize_tol
                if prt==1
                    println("t = $gam")
                end
                s = gam*p;z = z+s;x = z[1:n];y = z[n+1:n+m]
                F_old = F
                F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
                Y = F - F_old # Y is y in my notes
                append!(normFAll, normF)
                Zs = hcat(Zs,z)
                normF = LinearAlgebra.norm(F)
             else # null step
                if prt==1
                    println("null step")
                end
                s = p;z_temp = z+s;x_temp = z_temp[1:n];y_temp = z_temp[n+1:n+m]
                F_temp = [obj.∇xL(x_temp,y_temp); -obj.∇yL(x_temp,y_temp)]
                Y = F_temp - F_old
                append!(normFAll, normF)
                Zs = hcat(Zs,z)
            end
        else # no line search
            if prt==1
                println("fixed step = $fixed_stepsz")
            end
            gam=fixed_stepsz
            ###########################################
            # New edition on Jan 26, 2022, with fixed-stepsize method: use 1 after normF is tiny
            ###########################################
            if normF < 0.1
                #println("Broyden fixed-step switiching stepsize")
                gam = 1
            end
            s = gam*p;z = z+s;x = z[1:n];y = z[n+1:n+m]
            F_old = F
            F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
            Y = F - F_old # Y is y in my notes
            append!(normFAll, normF);Zs = hcat(Zs,z)
            normF = LinearAlgebra.norm(F)
        end
        stH = s'*H
        stHy = stH*Y
        sty = s'*Y
        H = H + ((s-H*Y)*stH)/stHy

        val = obj.L(x,y)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y))
        iter=iter+1
        if prt==1
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            println("||F|| = $normF")
            println("#################################End of iter $iter")
        end
    end#while
    append!(normFAll, normF)
    return x,y,iter,normFAll, val,normF,Zs
end

#######################################
#J-symm w and w/o line-search. The inverse Jacobian is obtained from  Sherman woodbury Morrison(SWM) to Hessian
#######################################
function secant_inv(x,y,obj,fixed_stepsz,sp, itNum,prt, use_ls,Ftol)
    val, ngx, ngy =0,0,0
    println("dummy")
    println("Secant method with using line-search = $use_ls")
    n = length(x);m = length(y);
    # J = [I(n)       zeros(n,m)
    #     zeros(m,n)     -I(m)]
    H = I(n+m);
    F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
    F_old = F
    normF=LinearAlgebra.norm(F)
    normFAll=[];z = [x;y];Zs=z
    iter=0
    while ( normF> Ftol ) && iter < itNum
        p = -H*F
        if(use_ls ==1)# line-search is being used
            gama = lin_search(x,y,obj,F,normF,p,c1)
            if gama >= stepsize_tol
                if prt==1
                    println("t = $gama")
                end
                s = gama*p;z = z+s;x = z[1:n];y = z[n+1:n+m]
                F_old = F;F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
                Y = F - F_old # Y is gradient differences
                append!(normFAll, normF);Zs = hcat(Zs,z);
                normF = LinearAlgebra.norm(F)
             else # null step
                if prt==1
                    println("null step")
                end
                gama=1
                F_old = F;s = p;z_temp = z+s;x_temp = z_temp[1:n];y_temp = z_temp[n+1:n+m]
                F_temp = [obj.∇xL(x_temp,y_temp); -obj.∇yL(x_temp,y_temp)]
                Y = F_temp - F_old
                append!(normFAll, normF);Zs = hcat(Zs,z)
            end
        else #line-search is not being used
            if prt==1
                println("fixed step = $fixed_stepsz")
            end
            gama=fixed_stepsz
            ###########################################
            # New edition on Jan 26, 2022, with fixed-stepsize method: use 1 after normF is tiny
            ###########################################
            if normF < 0.1
                #println("Jsymm fixed-step switiching stepsize")
                gama = 1
            end
            s = gama*p;z = z+s;x = z[1:n];y = z[n+1:n+m];
            F_old = F;F = [obj.∇xL(x,y); -obj.∇yL(x,y)]
            Y = F - F_old # Y is y in my notes
            append!(normFAll, normF);Zs = hcat(Zs,z)
            normF = LinearAlgebra.norm(F)
        end

        # apply J-symm update
        r = Y +gama*F_old # r = Y - D*s # i.e: Ds is just gamma*F_old
        ns2 = LinearAlgebra.norm(s)^2
        Js = [s[1:n];-s[n+1:n+m]];#J*s
        Jr = [r[1:n];-r[n+1:n+m]];#J*r
        α = s'*Jr
        a = r - α*Js/ns2
        Ha = H*a
        denom1 =ns2 + (s'*Ha)
        Ainv = H - (Ha*s'*H)/denom1
        AinvJs = Ainv*Js
        denom2 =ns2 + (Jr'*AinvJs)
        H = Ainv - (AinvJs*Jr'*Ainv)/denom2

        # print:
        val =obj.L(x,y)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y))
        iter=iter+1
        if prt==1
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            println("||F|| = $normF")
            println("#################################End of iter $iter")
        end
    end
    append!(normFAll, normF)

    return x,y,iter,normFAll, val, normF,Zs
end


##################################################################
# Trust-region method
##################################################################

function tr_dogleg(x,y, obj,sp, itNum,prt, F_tol, Del, max_Del, eta)
    n = length(x);m = length(y);val, ngx, ngy =0,0,0
    F = [obj.∇xL(x,y); -obj.∇yL(x,y)];normF=norm(F);
    nabla_F  = obj.∇F(x,y);
    B_k = I(n+m);H_k = I(n+m);

    g = nabla_F'*F;
    k=0    # count all iterations
    it=0   # count only actual steps
    normFAll=[];z = [x;y];Zs=z;

    while ( normF> F_tol ) && k < itNum
        pB = -H_k*(H_k'*g)                 #"compute quasi-Newton step"
        if( norm(pB) <= Del)
            s_k= pB
        else                               #compute the Cauchy point
            Bg = B_k*g
            norm_g = norm(g)
            Del2n_g=Del/norm_g
            ps = -Del2n_g*g
            norm_Bg_sqr = sum(Bg.*Bg)

            n_g2n_Bg=(norm_g)^2/norm_Bg_sqr
            if n_g2n_Bg >= Del2n_g        # Cauchy point outside
                s_k = ps
            else                          # Cauchy point inside, move along dogleg
                pU = -n_g2n_Bg*g
                q = pB-pU
                alpha = getAlpha(pU,q,Del)
                s_k = pU + alpha*q
            end
        end

        # Compue the prediction: ρ
        z_new = z+s_k
        F_new = [obj.∇xL(z_new[1:n],z_new[n+1:n+m]); -obj.∇yL(z_new[1:n],z_new[n+1:n+m])]
        numer = 0.5*(norm(F)^2 - norm(F_new)^2)
        Bs = B_k*s_k
        normBs2 = sum(Bs.*Bs)
        denom = -sum(g.*s_k) - 0.5*normBs2 # sum(g.*s_k)=g'*s_k
        if denom ==0
            rho = 1e99                 #division by 0
        else
            rho = numer/denom
        end

        # Expand or contract current trust-region radius: Δ_k
        if rho <= 0.5                 #contract
            Del = 0.5*Del
        else                          #expand
            Del = min(2*Del, max_Del)
        end

        Y = F_new - F

        # Update:
        r = Y - Bs
        norms2 = norm(s_k)^2
        Js = [s_k[1:n];-s_k[n+1:n+m]];#J*s_k
        Jr = [r[1:n];-r[n+1:n+m]];#J*r
        scalar0 = (Js'*r)/(norms2)^2;
        β = 1 # we prefer β=1 whenever B_{k+1} is nonsingular, which is almost surely the case
        update1= β*(Js*Jr'+ r*s_k')/norms2
        update2= β^2*(scalar0*Js*s_k')
        B_k = B_k + update1 - update2;
        H_k = B_k\I(m+n)

        # Check null step or valid step:
        if rho >=eta # sufficient decrease
            z = z_new;x = z_new[1:n];y = z_new[n+1:n+m];F = F_new
            normF = norm(F);nabla_F  = obj.∇F(x,y);g = nabla_F'*F
            val = norm(obj.L(x,y));ngx = norm(obj.∇xL(x,y));ngy = norm(obj.∇yL(x,y))
            it = it+1
        else
            if prt==1
                println("the step was null")
            end
        end
        append!(normFAll, normF)
        Zs = hcat(Zs,z)
        k=k+1
        if prt==1
            println("L = $val")
            println("|∇xL| = $ngx")
            println("|∇yL| = $ngy")
            println("||F|| = $normF")
            if(prob_type ==3)
                oval = obj.f(x)
                println("f(x) = $oval")
            end
            println("#################################End of iter $k")
        end
    end#while loop
    return x,y,k,normFAll, val, normF,Zs
end
