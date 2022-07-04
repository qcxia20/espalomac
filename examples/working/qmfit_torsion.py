import sys
sys.path[0] = "/pubhome/qcxia02/work/espalomac"
import espaloma as esp
import math
import torch
from tqdm import tqdm
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path, PosixPath
import argparse
from ase.units import kcal, mol, Hartree
if torch.cuda.is_available():
    from manager_torch import GPUManager
    gm = GPUManager()
    torch.cuda.set_device(gm.auto_choice(mode=0))
    device = torch.device("cuda")
else:
    device = torch.device("cpu")

kcalmol_per_Hartree =  Hartree / (kcal/mol)

def define_model():
    representation = esp.nn.Sequential( # GNN
        layer=esp.nn.layers.dgl_legacy.gn("SAGEConv"), # use SAGEConv implementation in DGL
        config=[128, "relu", 128, "relu", 128, "relu"], # 3 layers, 128 units, ReLU activation
    )
    readout = esp.nn.readout.janossy.JanossyPooling( # Pooling -> parameters
        in_features=128, config=[128, "relu", 128, "relu", 128, "relu"],
        out_features={              # define modular MM parameters Espaloma will assign
            4: {"k": 6}, # torsion barrier heights (can be positive or negative)
        },
    )
    espaloma_model = torch.nn.Sequential(
                    representation, readout,
                    esp.mm.geometry.GeometryInGraph(), # get torsion angle info
                    esp.mm.energy.EnergyInGraph(terms=["n4"]), # calculate energy
    )
    
    if torch.cuda.is_available():
        espaloma_model = espaloma_model.cuda()

    return espaloma_model

def train(bestpt,lastpt,lossfn,Nepoch,imgpath,args):
    # optimizer
    initial_lr = args.learningrate
    optimizer = torch.optim.Adam(espaloma_model.parameters(), initial_lr)
    # schedualer
    # scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.7,
                                                            # patience=5, min_lr=initial_lr/100) # copied from GeoMol

    train_losss, val_losss = [], []
    if (Path(modelpath) / bestpt).exists():
        pass # well-trained
    else:
        best_val_loss = math.inf
        best_epoch = 0
        # train
        for idx_epoch in tqdm(range(Nepoch)):
            losslist = []
            for g in ds_tr:
            # for g in ds:
                optimizer.zero_grad()
                g.heterograph = g.heterograph.to(device)
                g = espaloma_model(g.heterograph)
                loss = lossfn(g)
                loss.backward()
                optimizer.step()
                losslist.append(float(loss))
            
            train_losss.append(sum(losslist)/len(losslist))

        # val
            losslist = []
            with torch.no_grad():
                for g in ds_vl:
                    g.heterograph = g.heterograph.to(device)
                    g = espaloma_model(g.heterograph)
                    loss = loss_fn1(g)
                    losslist.append(float(loss))
            val_loss = sum(losslist)/len(losslist)
            val_losss.append(val_loss)

            if val_loss <= best_val_loss:
                best_val_loss = val_loss
                torch.save(espaloma_model.state_dict(), str(modelpath / bestpt))
            torch.save({
                'epoch': idx_epoch,
                'model': espaloma_model.state_dict(),
                'optimizer': optimizer.state_dict(),
                # 'scheduler': scheduler.state_dict(),
            }, str(modelpath / lastpt)) 


    train_losss = np.array(train_losss) # Ha -> kcal/mol 627.5
    val_losss = np.array(val_losss)

    ## Hartree to kcal/mol
    # train_losss = np.array(train_losss) * kcalmol_per_Hartree # Ha -> kcal/mol 627.5
    # val_losss = np.array(val_losss) * kcalmol_per_Hartree

    plt.plot(train_losss, label="train")
    plt.plot(val_losss, label="valid")
    plt.legend()
    plt.savefig(imgpath)
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Using truncated espaloma and charmm to fit torsion parameters")
    parser.add_argument("--Nepochs", type=int, help="number of epochs")
    parser.add_argument("--datasetpath", type=str, help="absolute path of dataset root path")
    parser.add_argument("--learningrate", type=float, default=1e-3, help="learning rate, e.g., 1e-3")
    args = parser.parse_args()

    N = args.Nepochs
    N1=N2 = N 
    # modelpath = Path("/pubhome/qcxia02/work/espalomac/examples/models") / f"{Path(args.datasetpath).name}_E{args.Nepochs}_L{args.learningrate}"
    # modelpath = Path("/pubhome/qcxia02/work/espalomac/examples/models") / f"single_E{args.Nepochs}_L{args.learningrate}"
    modelpath = Path("/pubhome/qcxia02/work/espalomac/examples/models") / f"single_E{args.Nepochs}_{Path(args.datasetpath).name}"
    if not modelpath.exists():
        modelpath.mkdir()
    img1path = modelpath / "train_step1.png"
    img2path = modelpath / "train_step2.png"

    imgpath = modelpath / "test.png"
    # def load_dataset():
    # dataset_name = "gen2"
    # dataset_name = "pepconf"
    # dataset_name = "vehicle"
    # dataset_name = "phalkethoh"
    dataset_name = args.datasetpath
    ds = esp.data.dataset.GraphDataset.load(dataset_name)
    ds.shuffle(seed=2022)
    ds_tr, ds_vl, ds_test = ds.split([8, 1, 1])
    ds_tr.save(str(modelpath/"train_dataset.pkl"))
    ds_vl.save(str(modelpath/"valid_dataset.pkl"))
    ds_test.save(str(modelpath/"test_dataset.pkl"))

    espaloma_model = define_model()
        
    # loss function
    loss_fn1 = esp.metrics.GraphMetric(
        # base_metric=torch.nn.MSELoss(), # use mean-squared error loss
        base_metric=torch.nn.MSELoss(), # use absolute mean-squared error loss, to make the energy value close to it
        between=['u', "u_charmm_totE"],         # between predicted and MM energies
        level="g", # compare on graph level
    )

    loss_fn2 = esp.metrics.GraphMetric(
        # base_metric=torch.nn.MSELoss(), # use mean-squared error loss
        base_metric=esp.metrics.center(torch.nn.MSELoss()), # use centered mean-squared error loss
        between=['u', "u_ref"],         # between predicted and QM energies
        level="g", # compare on graph level
    )

    # Train1
    train(bestpt="best1.pt",lastpt="last1.pt",lossfn=loss_fn1,Nepoch=N1,imgpath=img1path,args=args)

    # Train2
    train(bestpt="best2.pt",lastpt="last2.pt",lossfn=loss_fn2,Nepoch=N2,imgpath=img2path,args=args)

    # # 
    # # Inspection
    # inspect_metric = esp.metrics.center(torch.nn.L1Loss(), dim=0) # use mean-squared error loss
    # with torch.no_grad():
    #     # for idx_epoch in tqdm(range(N)):
    #     espaloma_model.load_state_dict(
    #         # torch.load("%s.th" % idx_epoch)
    #         # torch.load(str(modelpath / ("%s.th" % idx_epoch)), map_location=device)
    #         torch.load(str(modelpath / "best.pt"))
    #     )
    #     # training set performance
    #     u = []
    #     u_ref = []
    #     for g in ds_test:
    #         g.heterograph = g.heterograph.to(device)
    #         espaloma_model(g.heterograph)
    #         # u.append(g.nodes['g'].data['u'])
    #         u.append(g.nodes['g'].data['u'][0])
    #         # print(u)
    #         u_ref.append(g.nodes['g'].data['u_ref'][0])
    #         # print(u_ref)
    #     u = torch.cat(u, dim=0)
    #     u_ref = torch.cat(u_ref, dim=0)
    #     print("Best pt on test set:\n")

    #     print(float(inspect_metric(u, u_ref).to(torch.device("cpu"))))
    #     ## Hartree to kcal/mol
    #     # print(float(inspect_metric(u, u_ref).to(torch.device("cpu"))) * kcalmol_per_Hartree)
