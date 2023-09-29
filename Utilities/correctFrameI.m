function correctTimeIdx = correctFrameI(refArray, timeIdx1, numPoints)
corrArray = 1:0.5:floor(max(refArray)/2);
correctTimeIdx1 = floor(corrArray(timeIdx1));
correctTimeIdx2 = correctTimeIdx1+numPoints-1;
correctTimeIdx = correctTimeIdx1:correctTimeIdx2;
end