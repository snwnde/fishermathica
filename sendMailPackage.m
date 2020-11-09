(* ::Package:: *)

(* Original code by Rodrigo Murta, http://mathematica.stackexchange.com/questions/30836/encrypt-the-password-used-for-sendmail*)

BeginPackage["sendMailPackage`"]

sendMail::usage="sendMail[subject, body, to, file] send mail using myMail";


(*dir = DirectoryName[$InputFileName];
SetDirectory[dir];
*)


Begin["`Private`"]

pa=Uncompress@Import["myPass.txt","String"]

sendMail[subject_,body_,to_,file_:None]:=
SendMail[
    "To"->to,
    "Subject"->subject,
    "Body"->body,
    "From"->"casas@thphys.uni-heidelberg.de",
    "Server"->"extmail.urz.uni-heidelberg.de",
    "PortNumber"->587,
    "EncryptionProtocol"->"StartTLS",
    "UserName"->"ez378",
    "Password"->pa, (*Use compress to create it*)
    "Attachments"->file
]

SetAttributes[sendMail, {ReadProtected,Locked}];

End[]
EndPackage[]
