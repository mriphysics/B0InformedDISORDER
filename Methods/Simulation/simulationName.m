

function [fileName] = simulationName(rec)

fileName = '';

fileName = strcat(fileName, sprintf('SynthX(%d)T(%d)D(%d)_', ...
                            rec.Sim.synthX,...
                            rec.Sim.synthT,...
                            rec.Sim.synthD));
                        
fileName = strcat(fileName, sprintf('EstX(%d)T(%d)D(%d)_', ...
                            rec.Sim.estX,...
                            rec.Sim.estT,...
                            rec.Sim.estD));
                        
fileName = strcat(fileName, sprintf('ProvX(%d)T(%d)D(%d)_', ...
                            rec.Sim.provideX,...
                            rec.Sim.provideT,...
                            rec.Sim.provideD));
                        
fileName = strcat(fileName, sprintf('SNR(%d)_', rec.Sim.snrdB));    
fileName = strcat(fileName, sprintf('groupSw(%d-%d)_', rec.Sim.groupSweepsSynth,rec.Alg.parXT.groupSweeps));                        
fileName = strcat(fileName, sprintf('motion(%s)_', rec.Sim.motionType));                        
fileName = strcat(fileName, sprintf('NinterT(%d)_', rec.Sim.NinterleavesT));                        
fileName = strcat(fileName, sprintf('tran(%d)_', rec.Sim.tran));                        
fileName = strcat(fileName, sprintf('rot(%d)_', rec.Sim.rot));                        
fileName = strcat(fileName, sprintf('resol(%.2f-%.2f)', rec.Sim.resSynth,rec.Sim.resRecon(end)));
fileName = strcat(fileName, sprintf('.mat'));

