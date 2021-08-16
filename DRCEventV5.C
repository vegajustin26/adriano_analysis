#include "DRCEventV5.h"
#include <vector>


ClassImp(DRCHit)
ClassImp(DRCEvent)
ClassImp(DRCEventV5)

TClonesArray *DRCEvent::fgDRCHits = 0;

//______________________________________________________________________________
DRCEvent::DRCEvent()
{
  if (!fgDRCHits) fgDRCHits = new TClonesArray("DRCHit", 10);
  DRCHits = fgDRCHits;
}

//______________________________________________________________________________
DRCEvent::~DRCEvent()
{
  Reset();
}

//______________________________________________________________________________
void DRCEvent::Build(Int_t SpillNumber, Int_t TrigNum, Double_t eventTime, vector <Hit> drchits)
{
  //Build one event

  Reset();
  fSpillNumber   = SpillNumber;
  fTriggerNumber = TrigNum;
  fEventTime     = eventTime;
  DRCHit *drchit;
  fNDRCHits  = 0;

  for(UInt_t i1=0;i1<drchits.size();i1++)
  {
    drchit = AddDRCHit();
    drchit->fW.Set(drchits[i1].fW.size());
    drchit->fTimeSinceLastTrig=drchits[i1].fTimeSinceLastTrig;
    drchit->fCh=drchits[i1].fCh;
    for(UInt_t i2=0;i2<drchits[i1].fW.size();i2++)
    {
      drchit->fW.AddAt(drchits[i1].fW[i2],i2);
    }
  }
}

DRCHit *DRCEvent::AddDRCHit()
{
  TClonesArray &DRChits = *DRCHits;
  DRCHit *drchit = new(DRChits[fNDRCHits++]) DRCHit();
  return drchit;
}

void DRCEvent::Clear(Option_t *option)
{
  DRCHits->Clear(option);
}

void DRCEvent::Reset(Option_t *)
{
  // Static function to reset all static objects for this event

  DRCHits->Delete(/*option*/);
  fgDRCHits->Delete(/*option*/);
}




TClonesArray *DRCEventV5::fgDRCHitsV5 = 0;

//______________________________________________________________________________
DRCEventV5::DRCEventV5()
{
  if (!fgDRCHitsV5) fgDRCHitsV5 = new TClonesArray("DRCHit", 10);
  DRCHits = fgDRCHitsV5;
}

//______________________________________________________________________________
DRCEventV5::~DRCEventV5()
{
  Reset();
}

//______________________________________________________________________________
void DRCEventV5::Build(Int_t SpillNumber, Int_t MB1TrigNum, Int_t MB2TrigNum, Double_t eventTime, vector <Hit> drchits)
{
  //Build one event

  Reset();
  fSpillNumber      = SpillNumber;
  fMB1TriggerNumber = MB1TrigNum;
  fMB2TriggerNumber = MB2TrigNum;
  fEventTime        = eventTime;
  DRCHit *drchit;
  fNDRCHits  = 0;

  for(UInt_t i1=0;i1<drchits.size();i1++)
  {
    drchit = AddDRCHit();
    drchit->fW.Set(drchits[i1].fW.size());
    drchit->fTimeSinceLastTrig=drchits[i1].fTimeSinceLastTrig;
    drchit->fCh=drchits[i1].fCh;
    for(UInt_t i2=0;i2<drchits[i1].fW.size();i2++)
    {
      drchit->fW.AddAt(drchits[i1].fW[i2],i2);
    }
  }
}

DRCHit *DRCEventV5::AddDRCHit()
{
  TClonesArray &DRChits = *DRCHits;
  DRCHit *drchit = new(DRChits[fNDRCHits++]) DRCHit();
  return drchit;
}

void DRCEventV5::Clear(Option_t *option)
{
  DRCHits->Clear(option);
}

void DRCEventV5::Reset(Option_t *)
{
  // Static function to reset all static objects for this event

  DRCHits->Delete(/*option*/);
  fgDRCHitsV5->Delete(/*option*/);
}


