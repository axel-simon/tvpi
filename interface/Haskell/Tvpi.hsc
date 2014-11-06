-- Tvpi.hsc - Haskell binding for the TVPI domain
module Abs.Domain.Tvpi.Tvpi(
  Tvpi,
  Variable,
  new,
  join,
  copy,
  makeRoom,
  project,
  update,
  augment,
  getNoOfVariables,
  createUnboundedVariable,
  createVariable,
  Inequality(..),
  addInequalities,
  widen,
  entails,
  LinComponent(..),
  ApproximationResult(..),
  approximateInequality,
  dump
  ) where

import Foreign hiding (new)
import Foreign.C.String
import Monad (liftM)

#include "tvpi_c.h"

newtype Tvpi = Tvpi (ForeignPtr Tvpi)

type Variable = #type size_t

new :: Variable -> IO Tvpi
new vars = do
  tvpiPtr <- tvpiNew vars
  liftM Tvpi $ newForeignPtr tvpiPtr tvpiFree

foreign import ccall unsafe "tvpiNew" tvpiNew :: Variable -> IO (Ptr Tvpi)
foreign import ccall unsafe "&tvpiFree" tvpiFree :: FinalizerPtr Tvpi

join :: Tvpi -> Tvpi -> IO Tvpi
join (Tvpi tvpi1) (Tvpi tvpi2) =
  withForeignPtr tvpi1 $ \tvpiPtr1 ->
  withForeignPtr tvpi2 $ \tvpiPtr2 -> do
  tvpiPtr <- tvpiJoin tvpiPtr1 tvpiPtr2
  liftM Tvpi $ newForeignPtr tvpiPtr tvpiFree

foreign import ccall unsafe "tvpiJoin"
  tvpiJoin :: Ptr Tvpi -> Ptr Tvpi -> IO (Ptr Tvpi)

copy :: Tvpi -> IO Tvpi
copy (Tvpi tvpi) = withForeignPtr tvpi $ \tvpiPtr -> do
  tvpiPtr' <- tvpiCopy tvpiPtr
  liftM Tvpi $ newForeignPtr tvpiPtr' tvpiFree

foreign import ccall unsafe "tvpiCopy"
  tvpiCopy :: Ptr Tvpi -> IO (Ptr Tvpi)

makeRoom :: Tvpi -> Variable -> IO ()
makeRoom (Tvpi tvpi) vars = withForeignPtr tvpi $ \tvpiPtr -> do
  tvpiMakeRoom tvpiPtr vars

foreign import ccall unsafe "tvpiMakeRoom"
  tvpiMakeRoom :: Ptr Tvpi -> Variable -> IO ()

project :: Tvpi -> [Variable] -> Variable -> IO ()
project (Tvpi tvpi) vars headroom =
  withForeignPtr tvpi $ \tvpiPtr ->
  withArray vars $ \arrPtr ->
  tvpiProject tvpiPtr (fromIntegral $ length vars) arrPtr headroom

foreign import ccall unsafe "tvpiProject"
  tvpiProject :: Ptr Tvpi -> #{type size_t} -> Ptr Variable -> 
		  Variable -> IO ()

update :: Tvpi -> Variable -> Bool -> IO ()
update (Tvpi tvpi) var additive = withForeignPtr tvpi $ \tvpiPtr ->
  tvpiUpdate tvpiPtr var (fromBool additive)

foreign import ccall unsafe "tvpiUpdate"
  tvpiUpdate :: Ptr Tvpi -> Variable -> #{type int} -> IO ()

augment :: Tvpi -> Variable -> Variable -> IO ()
augment (Tvpi tvpi) target source = withForeignPtr tvpi $ \tvpiPtr ->
  tvpiAugment tvpiPtr target source

foreign import ccall unsafe "tvpiAugment"
  tvpiAugment :: Ptr Tvpi -> Variable -> Variable -> IO ()

getNoOfVariables :: Tvpi -> IO Variable
getNoOfVariables (Tvpi tvpi) = withForeignPtr tvpi $ tvpiGetNoOfVariables

foreign import ccall unsafe "tvpiGetNoOfVariables"
  tvpiGetNoOfVariables :: Ptr Tvpi -> IO Variable

createUnboundedVariable :: Tvpi -> IO Variable
createUnboundedVariable (Tvpi tvpi) = withForeignPtr tvpi $ \tvpiPtr ->
  tvpiCreateUnboundedVariable tvpiPtr

foreign import ccall unsafe "tvpiCreateUnboundedVariable"
  tvpiCreateUnboundedVariable :: Ptr Tvpi -> IO Variable

createVariable :: Tvpi -> Integer -> IO Variable
createVariable (Tvpi tvpi) val = withForeignPtr tvpi $ \tvpiPtr ->
  tvpiCreateVariable tvpiPtr (fromIntegral val)

foreign import ccall unsafe "tvpiCreateVariable"
  tvpiCreateVariable :: Ptr Tvpi -> #{type long} -> IO Variable

data Inequality = Inequality {
  a :: Integer,
  b :: Integer,
  c :: Integer }

addInequalities :: Tvpi -> Variable -> Variable -> [Inequality] -> IO Bool
addInequalities (Tvpi tvpi) varX varY ineqs =
  withForeignPtr tvpi $ \tvpiPtr -> do
    ineqPtrs <- mapM (\(Inequality a b c) -> inequalityNew (fromIntegral a)
			(fromIntegral b) (fromIntegral c)) ineqs
    withArray ineqPtrs $ \arrPtr ->
      tvpiAddInequalities tvpiPtr varX varY (fromIntegral $ length ineqPtrs)
        arrPtr 0

foreign import ccall unsafe "inequalityNew"
  inequalityNew :: #{type long} -> #{type long} -> #{type long} ->
		   IO (Ptr ())

foreign import ccall unsafe "tvpiAddInequalities"
  tvpiAddInequalities :: Ptr Tvpi -> Variable -> Variable ->
			 #{type size_t} -> Ptr (Ptr ()) -> #{type int} ->
			 IO Bool

widen :: Tvpi -> Tvpi -> IO ()
widen (Tvpi invar) (Tvpi changed) =
  withForeignPtr invar $ \invarPtr ->
  withForeignPtr changed $ \changedPtr ->
  tvpiWiden invarPtr changedPtr

foreign import ccall unsafe "tvpiWiden"
  tvpiWiden :: Ptr Tvpi -> Ptr Tvpi -> IO ()

entails :: Tvpi -> Tvpi -> IO Bool
entails (Tvpi tvpi1) (Tvpi tvpi2) =
  withForeignPtr tvpi1 $ \tvpi1Ptr ->
  withForeignPtr tvpi2 $ \tvpi2Ptr ->
  liftM toBool $ tvpiEntails tvpi1Ptr tvpi2Ptr

foreign import ccall unsafe "tvpiEntails"
  tvpiEntails :: Ptr Tvpi -> Ptr Tvpi -> IO #{type int}

newtype LinComponent = LinComponent (Integer, Variable)

instance Storable LinComponent where
  sizeOf _ = #{const sizeof(LinComponent)}
  alignment _ = alignment (undefined :: #{type long})
  peek ptr = do
    (c :: #{type long}) <- #{peek LinComponent, coefficient} ptr
    (v :: #{type Variable}) <- #{peek LinComponent, variable} ptr
    return (LinComponent (fromIntegral c, v))
  poke ptr (LinComponent (c, v)) = do
    #{poke LinComponent, coefficient} ptr (fromIntegral c :: #{type long})
    #{poke LinComponent, variable} ptr v

data ApproximationResult
  = NoInformation
  | Unsatisfiable
  | Satisfiable
  deriving Show

instance Enum ApproximationResult where
  toEnum (-1) = NoInformation
  toEnum 0 = Unsatisfiable
  toEnum 1 = Satisfiable
  fromEnum NoInformation = -1
  fromEnum Unsatisfiable = 0
  fromEnum Satisfiable = 1

approximateInequality :: Tvpi -> [LinComponent] -> Integer ->
			 IO ApproximationResult
approximateInequality (Tvpi tvpi) linComps const =
  withForeignPtr tvpi $ \tvpiPtr -> withArray linComps $ \arrayPtr ->
  liftM (toEnum . fromIntegral) $
    tvpiApproximateInequality tvpiPtr arrayPtr (fromIntegral (length linComps))
			      (fromIntegral const)

foreign import ccall unsafe "tvpiApproximateInequality"
  tvpiApproximateInequality :: Ptr Tvpi -> Ptr LinComponent ->
			       #{type size_t} ->
			       #{type long} -> IO (#{type int})

{-
dump :: Tvpi -> [String] -> IO ()
dump (Tvpi tvpi) varNames = withForeignPtr tvpi $ \tvpiPtr -> do
  strPtrs <- mapM newCString varNames
  withArray strPtrs $ tvpiPrint tvpiPtr (fromIntegral $ length strPtrs)
  mapM_ free strPtrs
-}

dump :: Tvpi -> IO ()
dump (Tvpi tvpi) = withForeignPtr tvpi tvpiShow

foreign import ccall unsafe "tvpiShow"
  tvpiShow :: Ptr Tvpi -> IO ()
