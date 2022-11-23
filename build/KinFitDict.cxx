// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME KinFitDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "/mnt/d/waleed/pp45/KinFit/KinFit/include/hdecaybuilder.h"
#include "/mnt/d/waleed/pp45/KinFit/KinFit/include/hkinfitter.h"
#include "/mnt/d/waleed/pp45/KinFit/KinFit/include/hneutralcandfinder.h"
#include "/mnt/d/waleed/pp45/KinFit/KinFit/include/hrefitcand.h"
#include "/mnt/d/waleed/pp45/KinFit/KinFit/include/hrootfitter.h"
#include "/mnt/d/waleed/pp45/KinFit/KinFit/include/hvertexfinder.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_HRefitCand(void *p = nullptr);
   static void *newArray_HRefitCand(Long_t size, void *p);
   static void delete_HRefitCand(void *p);
   static void deleteArray_HRefitCand(void *p);
   static void destruct_HRefitCand(void *p);
   static void streamer_HRefitCand(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HRefitCand*)
   {
      ::HRefitCand *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HRefitCand >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("HRefitCand", ::HRefitCand::Class_Version(), "hrefitcand.h", 13,
                  typeid(::HRefitCand), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HRefitCand::Dictionary, isa_proxy, 16,
                  sizeof(::HRefitCand) );
      instance.SetNew(&new_HRefitCand);
      instance.SetNewArray(&newArray_HRefitCand);
      instance.SetDelete(&delete_HRefitCand);
      instance.SetDeleteArray(&deleteArray_HRefitCand);
      instance.SetDestructor(&destruct_HRefitCand);
      instance.SetStreamerFunc(&streamer_HRefitCand);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HRefitCand*)
   {
      return GenerateInitInstanceLocal((::HRefitCand*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HRefitCand*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_HKinFitter(void *p);
   static void deleteArray_HKinFitter(void *p);
   static void destruct_HKinFitter(void *p);
   static void streamer_HKinFitter(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HKinFitter*)
   {
      ::HKinFitter *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::HKinFitter >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("HKinFitter", ::HKinFitter::Class_Version(), "hkinfitter.h", 53,
                  typeid(::HKinFitter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::HKinFitter::Dictionary, isa_proxy, 16,
                  sizeof(::HKinFitter) );
      instance.SetDelete(&delete_HKinFitter);
      instance.SetDeleteArray(&deleteArray_HKinFitter);
      instance.SetDestructor(&destruct_HKinFitter);
      instance.SetStreamerFunc(&streamer_HKinFitter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HKinFitter*)
   {
      return GenerateInitInstanceLocal((::HKinFitter*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HKinFitter*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *HNeutralCandFinder_Dictionary();
   static void HNeutralCandFinder_TClassManip(TClass*);
   static void delete_HNeutralCandFinder(void *p);
   static void deleteArray_HNeutralCandFinder(void *p);
   static void destruct_HNeutralCandFinder(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HNeutralCandFinder*)
   {
      ::HNeutralCandFinder *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::HNeutralCandFinder));
      static ::ROOT::TGenericClassInfo 
         instance("HNeutralCandFinder", "hneutralcandfinder.h", 23,
                  typeid(::HNeutralCandFinder), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &HNeutralCandFinder_Dictionary, isa_proxy, 0,
                  sizeof(::HNeutralCandFinder) );
      instance.SetDelete(&delete_HNeutralCandFinder);
      instance.SetDeleteArray(&deleteArray_HNeutralCandFinder);
      instance.SetDestructor(&destruct_HNeutralCandFinder);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HNeutralCandFinder*)
   {
      return GenerateInitInstanceLocal((::HNeutralCandFinder*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HNeutralCandFinder*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *HNeutralCandFinder_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::HNeutralCandFinder*)nullptr)->GetClass();
      HNeutralCandFinder_TClassManip(theClass);
   return theClass;
   }

   static void HNeutralCandFinder_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *HVertexFinder_Dictionary();
   static void HVertexFinder_TClassManip(TClass*);
   static void delete_HVertexFinder(void *p);
   static void deleteArray_HVertexFinder(void *p);
   static void destruct_HVertexFinder(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::HVertexFinder*)
   {
      ::HVertexFinder *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::HVertexFinder));
      static ::ROOT::TGenericClassInfo 
         instance("HVertexFinder", "hvertexfinder.h", 25,
                  typeid(::HVertexFinder), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &HVertexFinder_Dictionary, isa_proxy, 0,
                  sizeof(::HVertexFinder) );
      instance.SetDelete(&delete_HVertexFinder);
      instance.SetDeleteArray(&deleteArray_HVertexFinder);
      instance.SetDestructor(&destruct_HVertexFinder);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::HVertexFinder*)
   {
      return GenerateInitInstanceLocal((::HVertexFinder*)nullptr);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::HVertexFinder*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *HVertexFinder_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::HVertexFinder*)nullptr)->GetClass();
      HVertexFinder_TClassManip(theClass);
   return theClass;
   }

   static void HVertexFinder_TClassManip(TClass* ){
   }

} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr HRefitCand::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *HRefitCand::Class_Name()
{
   return "HRefitCand";
}

//______________________________________________________________________________
const char *HRefitCand::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HRefitCand*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int HRefitCand::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HRefitCand*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HRefitCand::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HRefitCand*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HRefitCand::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HRefitCand*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr HKinFitter::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *HKinFitter::Class_Name()
{
   return "HKinFitter";
}

//______________________________________________________________________________
const char *HKinFitter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HKinFitter*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int HKinFitter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::HKinFitter*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *HKinFitter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HKinFitter*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *HKinFitter::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::HKinFitter*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void HRefitCand::Streamer(TBuffer &R__b)
{
   // Stream an object of class HRefitCand.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TLorentzVector::Streamer(R__b);
      R__b >> cand;
      R__b >> fMomentum;
      R__b >> fTheta;
      R__b >> fPhi;
      R__b >> fR;
      R__b >> fZ;
      R__b >> fPid;
      fCov.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, HRefitCand::IsA());
   } else {
      R__c = R__b.WriteVersion(HRefitCand::IsA(), kTRUE);
      TLorentzVector::Streamer(R__b);
      R__b << cand;
      R__b << fMomentum;
      R__b << fTheta;
      R__b << fPhi;
      R__b << fR;
      R__b << fZ;
      R__b << fPid;
      fCov.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_HRefitCand(void *p) {
      return  p ? new(p) ::HRefitCand : new ::HRefitCand;
   }
   static void *newArray_HRefitCand(Long_t nElements, void *p) {
      return p ? new(p) ::HRefitCand[nElements] : new ::HRefitCand[nElements];
   }
   // Wrapper around operator delete
   static void delete_HRefitCand(void *p) {
      delete ((::HRefitCand*)p);
   }
   static void deleteArray_HRefitCand(void *p) {
      delete [] ((::HRefitCand*)p);
   }
   static void destruct_HRefitCand(void *p) {
      typedef ::HRefitCand current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_HRefitCand(TBuffer &buf, void *obj) {
      ((::HRefitCand*)obj)->::HRefitCand::Streamer(buf);
   }
} // end of namespace ROOT for class ::HRefitCand

//______________________________________________________________________________
void HKinFitter::Streamer(TBuffer &R__b)
{
   // Stream an object of class HKinFitter.

   TObject::Streamer(R__b);
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_HKinFitter(void *p) {
      delete ((::HKinFitter*)p);
   }
   static void deleteArray_HKinFitter(void *p) {
      delete [] ((::HKinFitter*)p);
   }
   static void destruct_HKinFitter(void *p) {
      typedef ::HKinFitter current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_HKinFitter(TBuffer &buf, void *obj) {
      ((::HKinFitter*)obj)->::HKinFitter::Streamer(buf);
   }
} // end of namespace ROOT for class ::HKinFitter

namespace ROOT {
   // Wrapper around operator delete
   static void delete_HNeutralCandFinder(void *p) {
      delete ((::HNeutralCandFinder*)p);
   }
   static void deleteArray_HNeutralCandFinder(void *p) {
      delete [] ((::HNeutralCandFinder*)p);
   }
   static void destruct_HNeutralCandFinder(void *p) {
      typedef ::HNeutralCandFinder current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HNeutralCandFinder

namespace ROOT {
   // Wrapper around operator delete
   static void delete_HVertexFinder(void *p) {
      delete ((::HVertexFinder*)p);
   }
   static void deleteArray_HVertexFinder(void *p) {
      delete [] ((::HVertexFinder*)p);
   }
   static void destruct_HVertexFinder(void *p) {
      typedef ::HVertexFinder current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::HVertexFinder

namespace {
  void TriggerDictionaryInitialization_libKinFitDict_Impl() {
    static const char* headers[] = {
"/mnt/d/waleed/pp45/KinFit/KinFit/include/hdecaybuilder.h",
"/mnt/d/waleed/pp45/KinFit/KinFit/include/hkinfitter.h",
"/mnt/d/waleed/pp45/KinFit/KinFit/include/hneutralcandfinder.h",
"/mnt/d/waleed/pp45/KinFit/KinFit/include/hrefitcand.h",
"/mnt/d/waleed/pp45/KinFit/KinFit/include/hrootfitter.h",
"/mnt/d/waleed/pp45/KinFit/KinFit/include/hvertexfinder.h",
nullptr
    };
    static const char* includePaths[] = {
"/mnt/d/waleed/pp45/KinFit/KinFit",
"/mnt/d/waleed/pp45/KinFit/KinFit/include",
"/home/waleed/miniconda3/envs/ROOT/include",
"/home/waleed/miniconda3/envs/ROOT/include/",
"/mnt/d/waleed/pp45/KinFit/KinFit/build/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libKinFitDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$hrefitcand.h")))  __attribute__((annotate("$clingAutoload$/mnt/d/waleed/pp45/KinFit/KinFit/include/hdecaybuilder.h")))  HRefitCand;
class __attribute__((annotate("$clingAutoload$hkinfitter.h")))  __attribute__((annotate("$clingAutoload$/mnt/d/waleed/pp45/KinFit/KinFit/include/hdecaybuilder.h")))  HKinFitter;
class __attribute__((annotate("$clingAutoload$/mnt/d/waleed/pp45/KinFit/KinFit/include/hneutralcandfinder.h")))  HNeutralCandFinder;
class __attribute__((annotate("$clingAutoload$/mnt/d/waleed/pp45/KinFit/KinFit/include/hvertexfinder.h")))  HVertexFinder;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libKinFitDict dictionary payload"

#ifndef __ROOFIT_NOBANNER
  #define __ROOFIT_NOBANNER 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "/mnt/d/waleed/pp45/KinFit/KinFit/include/hdecaybuilder.h"
#include "/mnt/d/waleed/pp45/KinFit/KinFit/include/hkinfitter.h"
#include "/mnt/d/waleed/pp45/KinFit/KinFit/include/hneutralcandfinder.h"
#include "/mnt/d/waleed/pp45/KinFit/KinFit/include/hrefitcand.h"
#include "/mnt/d/waleed/pp45/KinFit/KinFit/include/hrootfitter.h"
#include "/mnt/d/waleed/pp45/KinFit/KinFit/include/hvertexfinder.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"HKinFitter", payloadCode, "@",
"HNeutralCandFinder", payloadCode, "@",
"HRefitCand", payloadCode, "@",
"HVertexFinder", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libKinFitDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libKinFitDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libKinFitDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libKinFitDict() {
  TriggerDictionaryInitialization_libKinFitDict_Impl();
}
