function out = bitpack(data, bits, args)
    % data is logical ch x time. bits is ch x 1, out is 1 x T
    % assigns each channel of data into bits
    arguments
        data (:, :) logical
        bits (:, 1) uint16
        args.type = "uint16";
    end

    C = size(data, 1);
    assert(numel(bits) == C);

    G = zeros(1, C);
    for c = 1:C
        G(c) = 2^(bits(c) - 1);
    end

    out = cast(G * data, args.type); % 1 X t
    
end