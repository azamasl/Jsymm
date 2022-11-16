#######################################
# EGM with fixed stepsize for analytical center
######################################## 
function EGM(x,y,w,obj,fixed_stepsz,sp, itNum,prt,Ftol)
    BETA = 0.5
    println("EGM for analytic center with fixed stepsize")
    n = length(x);m = length(w);
    F = [obj.∇xL(x,y,w);
         obj.∇yL(x,y,w);
        -obj.∇wL(x,y,w)]
    F_old = F
    normF=LinearAlgebra.norm(F)
    normFAll=[normF];z = [x;y;w];
    Zs=z;iter=0

    val = obj.L(x,y,w)
    ngx = LinearAlgebra.norm(obj.∇xL(x,y,w))
    ngy = LinearAlgebra.norm(obj.∇yL(x,y,w))
    ngw = LinearAlgebra.norm(obj.∇wL(x,y,w))

    while ( normF> Ftol ) && iter < itNum
        p = -F
        gama = fixed_stepsz;
        while(minimum(y+gama*p[n+1:n+m]) <= 0)
            gama = BETA*gama;
        end
        Fk = [obj.∇xL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m]);
             obj.∇yL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m]);
            -obj.∇wL(x+gama*p[1:n], y+gama*p[n+1:n+m], w+gama*p[n+m+1:n+2*m])]

        ###updating z
        s = gama*p;zk = z+s;
        xk = zk[1:n];yk = zk[n+1:n+m];wk = zk[n+m+1:n+2*m];
        gama = fixed_stepsz;
        # Fk = [obj.∇xL(xk,yk,wk);
        #  obj.∇yL(xk,yk,wk);
        # -obj.∇wL(xk,yk,wk)]
        pk=-Fk
        while(minimum(y+gama*pk[n+1:n+m]) <= 0)
            gama = BETA*gama;
        end
        F = [obj.∇xL(x+gama*pk[1:n], y+gama*pk[n+1:n+m], w+gama*pk[n+m+1:n+2*m]);
             obj.∇yL(x+gama*pk[1:n], y+gama*pk[n+1:n+m], w+gama*pk[n+m+1:n+2*m]);
            -obj.∇wL(x+gama*pk[1:n], y+gama*pk[n+1:n+m], w+gama*pk[n+m+1:n+2*m])]
        s = gama*pk;z = z+s;
        x = z[1:n];y = z[n+1:n+m];w = z[n+m+1:n+2*m];
        
        normF = norm(F);
        #println(normF)
        append!(normFAll, normF);Zs = hcat(Zs,z)

        val = obj.L(x,y,w)
        ngx = LinearAlgebra.norm(obj.∇xL(x,y,w))
        ngy = LinearAlgebra.norm(obj.∇yL(x,y,w))
        ngw = LinearAlgebra.norm(obj.∇wL(x,y,w))
        iter=iter+1
    end
    append!(normFAll, normF)

    return x,y,w,iter,normFAll, val, normF,Zs
end


