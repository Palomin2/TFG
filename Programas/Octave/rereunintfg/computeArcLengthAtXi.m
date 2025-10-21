function s_val = computeArcLengthAtXi(Xi, p, knotVec, controlPts, weights, elRange, elConn) 
    noCtrPts = size(controlPts,1);
    noElems  = size(elConn,1);
    noGPsInt = max(6, p+3);      
    [Wint, Qint] = quadrature(noGPsInt, 'GAUSS', 1);
    
    s_val = 0;
    % bucle en elementos acumulando hasta eta
    for e = 1:noElems
        xiE   = elRange(e,:); % [xi_i,xi_i+1]
        conn  = elConn(e,:);
        wi=weights(conn);
        i = findspan(noCtrPts-1,p,Xi,knotVec);
        if Xi >= xiE(2)
            % [xi_e, xi_e1]
            for gp = 1:length(Wint)
                t = Qint(gp);
                xi_gp = 0.5*( (xiE(2)-xiE(1))*t + xiE(2) + xiE(1) );
                J2 = 0.5*(xiE(2)-xiE(1));
                N = BasisFuns(i+1, xi_gp, p, knotVec);      
                D = DerBasisFuns(i+1, xi_gp, p, 1, knotVec); 
                w = sum(N .* wi');
                W1 = sum(D(2,:) .* wi');
                dRdxi = ( (D(2,:) .* wi') * w - W1 * (N .* wi') ) / (w^2);
                jacob1 = dRdxi*controlPts(conn,1:2); 
                s_val = s_val + norm(jacob1) * J2 * Wint(gp);
            end
        else
            % entonces eta esta en este elemento:
            for gp = 1:length(Wint)
                t = Qint(gp);
                xi_gp = 0.5*( (Xi - xiE(1))*t + Xi + xiE(1) );  
                J2 = 0.5*(Xi - xiE(1)); 
                N = BasisFuns(i+1, xi_gp, p, knotVec);      
                D = DerBasisFuns(i+1, xi_gp, p, 1, knotVec); 
                w = sum(N .* wi');
                W1 = sum(D(2,:) .* wi');
                dRdxi = ( (D(2,:) .* wi') * w - W1 * (N .* wi') ) / (w^2);
                jacob1 = dRdxi*controlPts(conn,1:2); 
                s_val = s_val + norm(jacob1) * J2 * Wint(gp);
            end
            break; 
        end
    end
end