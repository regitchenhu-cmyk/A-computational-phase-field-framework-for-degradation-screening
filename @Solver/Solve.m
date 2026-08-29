function Solve(obj)
	%SOLVE solves a single time increment for the nonlinear system of
	%equations through a Newton-Raphson procedure, combined with a
	%staggered solution scheme

	outerStop = false;
	outerit = 0;
	err0_array = [];

	while outerStop==false		%Staggered scheme convergence loop
		outerit = outerit + 1;
		fprintf("Loop " + string(outerit) + ":\n");

		for stp=1:obj.NSteps	%staggered scheme steps
			fprintf("  SubStep " + string(stp) + "\n");
    		stop = false;
    		it = 0;
		
			stepnum = size(obj.convergence_log, 1)+1;
 		
			%perform once-per-step calculations
			obj.physics.OncePerStep(stp);

			%assemble system matrices
    		obj.physics.Assemble(stp);
    		obj.physics.Constrain(stp);


    		recalc_pre=true;
    		En_err0 = -1;
    		curr_max_it = obj.maxIt;
			while(stop == false)	%Newton-Raphson solver loop
		
        		fprintf("    Solving it:" + string(it) + "      ");
        		tsolve = tic;
        		
        		recalc_pre = true;
        		if (recalc_pre)
            		%[P,R,C] = equilibrate(obj.physics.K{stp});
            		%recalc_pre = false;
        		end
		
        		if false   %use preconditioned system
            		d = -R*P*obj.physics.fint{stp};
            		B = R*P*obj.physics.K{stp}*C;
					%cond_num(it+1)=condest(B);
					if true
						dy = B\d;
					else
						[L,U] = ilu(B,struct('type','nofill'));
						dy = gmres(B,d,[],1e-4,500,L,U);
					end
            		dx = C*dy;
				else %do not allow any preconditioning
            		dx = -obj.physics.K{stp}\obj.physics.fint{stp};
        		end
        		tsolve = toc(tsolve);
        		fprintf("        (Solver time:"+string(tsolve)+")\n");
				%fprintf("Conditioning numbers: "+string(cond_num(it+1))+"\n");
		
				%line-search
        		if (obj.linesearch && it>-1)
            		e0 = obj.physics.fint{stp}'*dx;
            		obj.physics.Update(dx, stp);
            		
            		obj.physics.Assemble(stp);
            		obj.physics.Constrain(stp);
            		
            		e1 = obj.physics.fint{stp}'*dx;
            		factor = -e0/(e1-e0);
            		factor = max(obj.linesearchLims(1), min(obj.linesearchLims(2), factor));
            		obj.physics.Update(-(1-factor)*dx, stp);
            		fprintf("    Linesearch: " + string(e0) + " -> " + string(e1) + ":  eta=" + string(factor) +"\n");
        		else
            		obj.physics.Update(dx, stp);
        		end
        		
        		% re-assemble system
        		obj.physics.Assemble(stp);
        		obj.physics.Constrain(stp);

				% convergence check
        		if (En_err0 < 0)
            		En_err0 = sum(abs(obj.physics.fint{stp}.*dx));
            		En_err = En_err0;
		
            		if (En_err0==0)
                		En_err0 = 1e-12;
					end

					err0_array(stp) = En_err0;
        		else
            		En_err = sum(abs(obj.physics.fint{stp}.*dx));
        		end
        		        		En_err_n = En_err/En_err0;
		
				if isnan(En_err_n) || isinf(En_err_n)
					error('Solver:NonFiniteResidual', ...
						'Residual is NaN or Inf in staggered substep %d; state was not committed.', stp);
				end
		
				obj.convergence_log(stepnum,stp,it+1) = En_err_n;
                    		
        		fprintf('    Residual: %e   (%e)\n', En_err_n, En_err);
        		
        		it=it+1;
				if (it>=curr_max_it || En_err_n<obj.Conv || En_err<obj.tiny)
					stop = true;
				end
			end

			if ~(En_err_n < obj.Conv || En_err < obj.tiny)
				error('Solver:NonConvergence', ...
					['Staggered substep %d did not converge within %d Newton iterations ' ...
					 '(relative residual %.3e, absolute residual %.3e); state was not committed.'], ...
					stp, curr_max_it, En_err_n, En_err);
			end
		end

		% Code-level stopping heuristic: the stored first residual-increment
		% products retain the native scaling and units of their subproblems.
		outer_err = max(abs(err0_array));
		if outer_err < obj.tiny
			outerStop = true;
		elseif outerit >= obj.OuterLoops
			error('Solver:OuterNonConvergence', ...
				['Staggered iteration did not converge within %d outer sweeps ' ...
				 '(outer residual %.3e); state was not committed.'], ...
				obj.OuterLoops, outer_err);
		end

	end
    
	% Commit only after all nonlinear and staggered convergence checks pass.
	obj.physics.Commit("Pathdep");
    obj.physics.Commit("Timedep");
end

