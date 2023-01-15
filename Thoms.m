function ps = Thoms(N,p,w,deltax,deltay)
ps = p;
for j = 1:1:N
        for i = 1:1:N
            if i == 2
                ps(i+(j-1)*N,1)=deltax^2*w(i+(j-1)*N-1,1)/2;
            else
                if i == N-1
                    ps(i+(j-1)*N,1)=deltax^2*w(i+(j-1)*N+1,1)/2;
                end
            end
            if j == 2
                ps(i+(j-1)*N,1)=deltay^2*w(i+(j-2)*N,1)/2;
            else if j == N-1
                    ps(i+(j-1)*N,1)=deltay^2*w(i+(j)*N,1)/2;
                end
            end
        end
    end
end