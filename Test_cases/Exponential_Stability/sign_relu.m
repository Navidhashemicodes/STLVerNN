function y=sign_relu(x,c)
w=1/c;
y=poslin(w*x+1)-poslin(w*x-1)-1;
end

    