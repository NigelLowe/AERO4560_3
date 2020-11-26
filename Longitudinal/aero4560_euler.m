function  [X_out] = aero4560_euler(DT, A, B, X, U)

    X1_dot = A * X + B * U; 
    Am = X1_dot * DT;                   
    X2_dot = A * (X+Am/2) + B * U; 
    Bm = X2_dot * DT;
    X3_dot = A * (X+Bm/2) + B * U;
    Cm = X3_dot * DT;
    X4_dot = A * (X+Cm/2) + B * U;
    Dm = X4_dot * DT;           
    X_out = X + (Am+2*Bm+2*Cm+Dm)/6;  
    
end



