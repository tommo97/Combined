function x = test_sort(x)

for i = 1:length(x)
    xi = x(i);
    for j = i+1:length(x)
        if (x(j)<xi)
            x(i) = x(j);
            x(j) = xi;
        end
    end
end