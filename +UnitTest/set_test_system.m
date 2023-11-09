function scl = set_test_system()
    dim = randi([20 100], 1);
    d = rand(2, 1) * 5 + 1;
    r = randi([0 4], 1);
    
    if mod(r, 4) == 0
        domain = [-d(1) d(1); -d(2) d(2)];
    elseif mod(r, 4) == 1
        domain = [-d(1) d(1);0 d(2)];
    elseif mod(r, 4) == 2
        domain = [0 d(1);-d(2) d(2)];
    else
        domain = [0 d(1);0 d(2)];
    end
            
    scl = SystemInfo(dim, domain);
end
    
