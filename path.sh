
case "$USER" in
"ac1nmx")
  # Data root
  #export DATA_ROOT="/fastdata/ac1nmx/sound-localisation"
  export DATA_ROOT="/data/ac1nmx/data/sound-localisation"
  ;;
"ning")
  # Data root
  export DATA_ROOT="/Users/ning/work/TwoEars/data/sound-localisation"
  ;;
*)
  echo "Please define WAV_ROOT and REC_ROOT for user $USER"
  ;;
esac

export PATH=$PWD/local:$PATH

