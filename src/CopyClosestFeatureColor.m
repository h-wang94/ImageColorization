function [lab_out, tiesIdx] = CopyClosestFeatureColor(source, samples, target, GRAPH)
%Color transfer technique that assigns for each target pixel the color the
%source pixel with smallest distance in feature space.

nCandidates = 10;
colors = samples.ab;

lab_out = zeros(size(target.image));
lab_out(:,:,1) = target.luminance*100;

[r, c] = size(target.luminance);
tiesIdx = [];
for j = 1:c
    for i = 1:r
        idx = (i-1) + (j-1)*r + 1;
        [best_idxs, ~, ties] = BestMatchesFS(target.fv(:, idx), target.fv_w, samples.fv, nCandidates);

        %1NN color transfer.
        lab_out(i, j, 2:3) = colors(:, best_idxs(1));
        %TODO: implementar media geometrica para chamar aqui

        %% Debug:
        if (ties > 1)
            tiesIdx = [tiesIdx [i j]'];
        end

        CANDIDATES = false;
        if (ties > 1 && CANDIDATES)
            subplot(1,2,1)
            imshow(lab2rgb(lab_out)); hold on;
            scatter(j, i, '*r'); hold off;
            subplot(1,2,2)
            imshow(source.image); hold on;
            img_idxs = samples.idxs(:,best_idxs);
            scatter(img_idxs(2,1), img_idxs(1,1), '*g');
            scatter(img_idxs(2,2:end), img_idxs(1,2:end), '*r'); hold off;
            drawnow
        end
    end
end


end

