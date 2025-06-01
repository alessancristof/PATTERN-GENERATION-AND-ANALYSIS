function varargout = structureDescriptors(FirstVariation,SecondVariation,FirstEigenvector1,FirstEigenvector2,varoutlist)
% Retrieve queried variables from structure tensor.

% INPUTS
% FirstVariation: square root of maximal eigenvalue of structure tensor.
% SecondVariation: square root of minimal eigenvalue of structure tensor.
% Vmax1: first component of eigenvector associated with maximal eigenvalue.
% Vmax2: second component of eigenvector associated with maximal eigenvalue.

L1=[];
L2=[];
Ltot=[];
Ldiff=[];

for i=1:numel(varoutlist)
    switch varoutlist{i}
        case 'First Variation'
            varargout{i}=FirstVariation;
        case 'Second Variation'
            varargout{i}=SecondVariation;
        case 'Total Variation' % square root of the trace
            if isempty(Ltot)
                if isempty(L1)
                    L1=FirstVariation.^2;
                end
                if isempty(L2)
                    L2=SecondVariation.^2;
                end
                Ltot=L1+L2;
            end
            varargout{i}=sqrt(Ltot);
        case 'Variation Difference'
            varargout{i}=FirstVariation-SecondVariation;
        case 'First Eigenvalue'
            if isempty(L1)
                L1=FirstVariation.^2;
            end
            varargout{i}=L1;
        case 'Second Eigenvalue'
            if isempty(L2)
                L2=SecondVariation.^2;
            end
            varargout{i}=L2;
        case 'Trace'
            if isempty(Ltot)
                if isempty(L1)
                    L1=FirstVariation.^2;
                end
                if isempty(L2)
                    L2=SecondVariation.^2;
                end
                Ltot=L1+L2;
            end
            varargout{i}=Ltot;
        case 'Linear Eccentricity'
            if isempty(Ldiff)
                if isempty(L1)
                    L1=FirstVariation.^2;
                end
                if isempty(L2)
                    L2=SecondVariation.^2;
                end
                Ldiff=max(L1-L2,0);
            end
            varargout{i}=sqrt(Ldiff);
        case 'Eccentricity'
            if isempty(L1)
                L1=FirstVariation.^2;
            end
            if isempty(Ldiff)
                if isempty(L2)
                    L2=SecondVariation.^2;
                end
                Ldiff=max(L1-L2,0);
            end
            varargout{i}=sqrt(Ldiff./(L1+eps));
        case 'Coherence'
            if isempty(Ltot)
                if isempty(L1)
                    L1=FirstVariation.^2;
                end
                if isempty(L2)
                    L2=SecondVariation.^2;
                end
                Ltot=L1+L2;
            end
            if isempty(Ldiff)
                if isempty(L1)
                    L1=FirstVariation.^2;
                end
                if isempty(L2)
                    L2=SecondVariation.^2;
                end
                Ldiff=max(L1-L2,0);
            end
            varargout{i}=(Ldiff./(Ltot+eps)).^2;
        case 'DA'
            varargout{i}=max(FirstVariation-SecondVariation,0)./(FirstVariation+eps);
        case 'First Angle' % in image coordinates reference system (not intrinsic)
            % note that components are swapped and changed sign to have 0 angle along x axis of image
            var=FirstEigenvector2;
            var(abs(var)<=eps)=eps;
            varargout{i}=atan(-FirstEigenvector1./var)/pi*180;
        case 'Second Angle' % in image coordinate reference system (not intrinsic)
            % note that components are swapped and changed sign to have 0 angle along x axis of image
            var=FirstEigenvector1;
            var(abs(var)<=eps)=eps;
            varargout{i}=atan(FirstEigenvector2./var)/pi*180;
        case 'First Eigenvector'
            varargout{i}=cat(4,FirstEigenvector1,FirstEigenvector2);
        case 'Second Eigenvector'
            varargout{i}=cat(4,-FirstEigenvector2,FirstEigenvector1);
    end
end



