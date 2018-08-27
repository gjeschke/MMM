function r=vec_dist(a1,e1,a2,e2)
%
% Distance between two vectors with start points a1, a2 and end points e1,
% e2
%
% (c) G. Jeschke, 2007

% test, if identical
if sum(abs(a1-a2))<10*eps && sum(abs(e1-e2))<10*eps,
    r=0;
else
    % points with shortest distance on unbounded straight lines
    % see script "Computer graphics", A. Kolb, Uni Siegen, slide 175 ff.
    % our a1, a2 are Kolb's V1,V2
    % direction vectors
    I1=e1-a1;
    I2=e2-a2;
    b11=dot(I1,I1);
    b22=dot(I2,I2);
    sc=dot(I1,I2);
    bflip=[b22, sc; sc, b11]; % the matrix used in calculating the alpha
    v=[dot(a2-a1,I1);dot(a1-a2,I2)];
    deter=bflip(1,1)*bflip(2,2)-sc^2;
    if deter~=0,
        nf=1/deter;
        alph=nf*bflip*v;
        % the points of shortest distance between the unbounded lines
        p1=a1+alph(1)*I1;
        p2=a2+alph(2)*I2;
        if alph(1)>=0 && alph(1)<=1,
            if alph(2)>=0 && alph(2)<=1,
                r=norm(p2-p1);
            else
                % clip alpha2
                if alph(2)<0, p2c=a2; else p2c=e2; end;
                a1c=dot(p2c-a1,I1)/b11;
                if a1c<0, a1c=0; end;
                if a1c>1, a1c=1; end;
                p1c=a1+a1c*I1;
                r=norm(p2c-p1c);
            end;
        elseif alph(2)>=0 && alph(2)<=1,
            % clip alpha2
            if alph(1)<0, p1c=a1; else p1c=e1; end;
            a2c=dot(p1c-a2,I2)/b22;
            if a2c<0, a2c=0; end;
            if a2c>1, a2c=1; end;
            p2c=a2+a2c*I2;
            r=norm(p2c-p1c);
        else
            % compute both clippings
            if alph(2)<0, p2c=a2; else p2c=e2; end;
            a1c=dot(p2c-a1,I1)/b11;
            if a1c<0, a1c=0; end;
            if a1c>1, a1c=1; end;
            p1c=a1+a1c*I1;
            r1=norm(p2c-p1c);
            if alph(1)<0, p1c=a1; else p1c=e1; end;
            a2c=dot(p1c-a2,I2)/b22;
            if a2c<0, a2c=0; end;
            if a2c>1, a2c=1; end;
            p2c=a2+a2c*I2;
            r2=norm(p2c-p1c);
            r=min(r1,r2);
        end;
    else
        % parallel lines
        lambda1=dot(a1-a2,I2)/norm(I2)^2;
        if lambda1<0, lambda1=0; end;
        if lambda1>1, lambda1=1; end;
        c1=a1-(a2+lambda1*I2);
        r1=norm(c1);
        lambda2=dot(e1-a2,I2)/norm(I2)^2;
        if lambda2<0, lambda2=0; end;
        if lambda2>1, lambda2=1; end;    
        c2=e1-(a2+lambda2*I2);
        r2=norm(c2);
        r=min(r1,r2);
    end;
end;