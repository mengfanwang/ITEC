function bndwrite(data, label, path)
    outIm = display_bnd_map_3d(data,label);
    write4dTiffRGB(uint8(outIm),[path '.tif']);
end