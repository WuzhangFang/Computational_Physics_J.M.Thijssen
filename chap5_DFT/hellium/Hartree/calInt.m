function result = calInt(arr, N, h)
    help = (arr(1)*17+arr(2)*59+arr(3)*43+arr(4)*49)/48;
    for i=5:N-4
        help = help + arr(i);
    end
    help = (arr(N)*17+arr(N-1)*59+arr(N-2)*43+arr(N-3)*49)/48+help;
    result = help * h;
end