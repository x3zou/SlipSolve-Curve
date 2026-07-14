function not_implemented(stageName, contract)
%NOT_IMPLEMENTED Shared placeholder for stage wrappers.

message = sprintf("%s is not implemented yet.\n\nExpected contract:\n%s", stageName, contract);
error("SlipSolve:StageNotImplemented", "%s", message);

end

