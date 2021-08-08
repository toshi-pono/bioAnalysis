window.onload = () => {
  const canvas = document.getElementById("canvas");
  const context = canvas.getContext("2d");

  context.fillRect(0, 0, 10, 10);

  const aminos = [
    "MKPIIAVNYKAYYPYSFGENALRIARDAKRVWEETGVEVILAPPFTEIYRVLKEVEGSGVKVFAQHAVMKDMALNALKAVS",
    "MAPRKFFVGGNWKMNGDKKSLGELIHTLNGAKLSADTEVVCGAPSIYLDFARQKLDAKIGVAAQNCYCKELASQHDVDGFLVGGASLKPEFVDIINAKH",
    "MRRYLIAGNWKMNTSLETGTALASGLADHVRGRDLPVDVLVCPPFPYLAAVKATAGEAGISVGAQNCYFEASGAFTGEVSVDMLKDIGCDS",
  ];
  const data = {
    aminos: aminos,
  };
  const getTree = async () => {
    console.log(JSON.stringify(data));
    const res = await fetch("/tree", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify(data),
    });
    const tree = await res.json();
    console.log(tree);
  };
  getTree();
};
