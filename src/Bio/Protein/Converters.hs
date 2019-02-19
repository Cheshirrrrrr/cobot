{-# OPTIONS_GHC -fno-warn-orphans  #-}
{-# LANGUAGE RecordWildCards  #-}
{-# LANGUAGE TypeApplications #-}

module Bio.Protein.Converters where

import           Bio.Molecule          (Molecule (..))
import           Bio.Protein.AminoAcid (AminoAcid (..), BBCAT, BBT, Env (..),
                                        FromThreeSymbols (..))
import           Bio.Protein.Chain     (ProteinChain (..))
import           Bio.Structure         (Atom (..), Chain (..), Model (..),
                                        Residue (..),
                                        StructureSerializable (..))
import           Bio.Utils.Geometry    (V3R)
import           Bio.Utils.IUPAC       (AtomType (..))
import           Control.Lens          (Const (..))
import           Data.Array            (Array, assocs, elems)
import           Data.Functor.Identity (Identity (..))
import           Data.List             (find)
import           Data.Maybe            (fromMaybe)
import           Data.Text             (unpack)

instance FromResidue a => StructureSerializable [Molecule Int (ProteinChain Int a)] where
    serializeModels = fmap (Molecule . assocs . fmap (ProteinChain . fmap fromResidue . chainResidues) . modelChains) . elems

findAtom :: AtomType -> [Atom] -> Atom
findAtom atomType = fromMaybe (error $ "No " <> show atomType <> " atom in model")
                  . find ((== atomType) . read @AtomType . unpack . atomName)

getCoords :: Array Int Atom -> AtomType -> Identity V3R
getCoords atoms = pure . atomCoords . flip findAtom (elems atoms)

class FromResidue a where
    fromResidue :: Residue -> a

instance FromResidue (BBT V3R) where
    fromResidue Residue{..} = AminoAcid (pure nAtom) (Env caAtom (Const . fromThreeSymbols $ resName)) (pure cAtom)
      where
        nAtom  = getCoords resAtoms N
        caAtom = getCoords resAtoms CA
        cAtom  = getCoords resAtoms C

instance FromResidue (BBCAT V3R) where
    fromResidue Residue{..} = AminoAcid (Const ()) (Env caAtom (Const . fromThreeSymbols $ resName)) (Const ())
      where
        caAtom = getCoords resAtoms CA
