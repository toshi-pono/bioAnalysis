window.onload = () => {
  const canvas = document.getElementById("canvas");
  const context = canvas.getContext("2d");

  // todo get by canvas
  const WIDTH = 640;
  const HEIGHT = 480;
  const LINE_WIDTH = 2;
  const PADDING = 100;

  const aminos = [
    "MKPIIAVNYKAYYPYSFGENALRIARDAKRVWEETGVEVILAPPFTEIYRVLKEVEGSGVKVFAQHADPVEPGAVTGYIPVEGLKEAGVHGVILNHSEHRLKIADINALIIKARRLGLKTLACADVPETGAAIALLKPDMIAVEPPELIGTGVSVSKAKPEVITNSVSMIRSVNKEALILTGAGITTGEDVYQAVKLGTIGVLVASGIVKAKDPYSVMKDMALNALKAVS",
    "MAPRKFFVGGNWKMNGDKKSLGELIHTLNGAKLSADTEVVCGAPSIYLDFARQKLDAKIGVAAQNCYKVPKGAFTGEISPAMIKDIGAAWVILGHSERRHVFGESDELIGQKVAHALAEGLGVIACIGEKLDEREAGITEKVVFEQTKAIADNVKDWSKVVLAYEPVWAIGTGKTATPQQAQEVHEKLRGWLKSHVSDAVAQSTRIIYGGSVTGGNCKELASQHDVDGFLVGGASLKPEFVDIINAKH",
    "MARKFFVGGNWKCNGTAEEVKKIVNTLNEAQVPSQDVVEVVVSPPYVFLPLVKSTLRSDFFVAAQNCWVKKGGAFTGEVSAEMLVNLDIPWVILGHSERRAILNESSEFVGDKVAYALAQGLKVIACVGETLEEREAGSTMDVVAAQTKAIADRVTNWSNVVIAYEPVWAIGTGKVASPAQAQEVHDELRKWLAKNVSADVAATTRIIYGGSVNGGNCKELGGQADVDGFLVGGASLKPEFIDIIKAAEVKKSA",
    "MRQIIIAGNWKMHKTQTESLEFLQGFLSHLEDTPEERETVLCVPFTCLNFMSKNLHGSRVKLGAQNIHWADQGAFTGEISGEMLKEFGINYVIVGHSERRQYFGETDETVNARLLAAQKHGLTPILCVGESKAQRDAGETEAVISAQIEKDLVNVDQNNLVIAYEPIWAIGTGDTCEAAEANRVIGLIRSQLTNKNVTIQYGGSVNPKNVDEIMAQPEIDGALVGGASLDPESFARLVNYQ",
    "MHKTQAESLEFLQSFLPQLENTAEDREVILCAPYTALGVMSKNLHGTRVRIGSQNVHWEESGAFTGEIAPSMLTEIGVTYAVVGHSERRQYFGETDETVNFRARAAQKAELTPILCVGESKEQRDAGQTETVIKEQLKADLVGVDLSQLVIAYEPIWAIGTGDTCEAEEANRVIGMIRSELSSSDVPIQYGGSVKPANIDEIMAQPEIDGALVGGASLDPVGFARIVNYEAT",
    "MHKTQAESLEFLQSFLPQLENTAEDREVILCAPYTALGVMSKNLHGTRVRIGSQNVHWEESGAFTGEIAPSMLTEIGVTYAVVGHSERRQYFGETDETVNFRARAAQKAELTPILCVGESKEQRDAGQTETVIKEQLKADLVGVDLSQLVIAYEPIWAIGTGDTCEAEEANRVIGMIRSELSSSDVPIQYGGSVKPANIDEIMAQPEIDGALVGGASLDPVGFARIVNYEAT",
    "MKRQIVIAGNWKMHKTNSEAMQLANQVRIKTMDITKTQIVICPPFTALAPVYEVIGDSRIHLGAQNMFWEKEGAFTGEISAGMIKSTGADYVIIGHSERRQYFGESDETVNKKVKAALENGLKPIVCVGETLEEREANITLKVVSRQIRGAFADLSAEQMKKVIVAYEPVWAIGTGKTATPEQAQQVHQEIRQLLTEMFGSEIGEKMVIQYGGSVKPANAESLLSQPDIDGALVGGACLKADSFSEIIHIAEKLQ",
    "MKTRQQIVAGNWKMNKNYGEGRELAMEIVERLKPSNTQVVLCAPYIHLQLVKNIIKDVASLYLGAQNCHQEDKGAYTGEISVDMLKSVGVSYVILGHSERREYFGESDELLAKKTDKVLAAGLLPIFCCGESLDIRDAGTHVAHVQAQIKAGLFHLSPEEFQKVVIAYEPIWAIGTGRTASPEQAQDMHAAIRALLTDQYGAEIADATTILYGGSVNGGNAAVLFSQPDVDGGLVGGASLKAEEFITIVEATKK",
    "MVYWVGTSWKMNKTLAEAMDFAAILAGFVPGFDDRIQPFVIPPFTAVRQVKQALSSTRVKVGAQNMHWADAGAWTGEISPVMLTDCGLDLVELGHSERREHFGETDRTVGLKTAAAVKHGLIPLICVGETLAERESGEADAVLAKQVEGALQFFEEEVKGATILFAYEPVWAIGDKGIPASSDYADKQQGLIKAVAGSLLPSVPSVLYGGSVNPGNAAELIGQPNVDGLFIGRSAWQAQGYIDILGRASAAI",
    "MRRYLIAGNWKMNTSLETGTALASGLADHVRGRDLPVDVLVCPPFPYLAAVKATAGEAGISVGAQNCYFEASGAFTGEVSVDMLKDIGCDSVILGHSERRHVIKECDDMINKKTKAAIEGGLQVVLCVGELLEEREADKTEAVLDEQMAGGLKDISAEQMTNVVIAYEPVWAIGTGKTASPEQAEQAHAHLRKWLADRYTSEVAEQTRILYGGSVKPANAKELLGQQNVDGALVGGASLTVDNFGPIIDAGVELSA",
    "MPEEKPVIMINFKTYNESYGFRAHDIAEAAETVAEESGIEIVICPGFMDIHPMSNHYRLPVFAQHIDGISPGAHTGHILAEAVRAAGATGTLINHSERRLTLADISAAVDAAKRANLKTVVCTNNTATSGAAAALSPDYVAIEPPELIGSGISVATADPEIIENSVNAVKSVNKDVKVLAGAGISSGSCVKRAVELGSDGVLLASGVVKAEDPAVVLRDLVSKI",
    "MGSPLIVVNFKTYLEGTGERSVDIARACRDVAEDSGVDIAVAPQMCDIYRVASMVDIPVYSQHVDGIGAGSFTGHAFAPAIKEAGASGTLINHSENRLTLADIEAAIQASKAVGLKTIVCTNNIPTSAAAAALSPDYVAVEPPELIGSGIPVSEADPDVVKGSVEAVMNIDSGVSVLCGAGISKGKDLKAALDLGSKGVLLASGIVKSEDPRSAMEDLISLI",
    "MRKKIVAGNWKMNLDYTEGLTLFSEVINMIKDEVTGSQQAVVCSPFIHLHSLVQLGKDYNKVSVGAQNAHQAEAGAYTGEISSRMIKSVGAEYVIIGHSERRQYFGETNDLLAKKTDAVLKNQLTPIFCIGETLQERETEKHFEVIKSQLLEGVFHLDETAFAKLVIAYEPVWAIGTGVTASAEQAQEIHAFIRAEIAQKYSQQVADDITILYGGSCNPKNAAELFAKGDIDGGLIGGASLKSRDFVDILKVFN",
    "MATSKTVGRVPLMAGNWKMNLDHLQATHLIQKLDWTLRDAKHDYDGVEVAVLPPFTDLRSVQTLVEGDRLHLRYGAQDLSPHASGAYTGDISGAFLKKLGCTYVVVGHSERREGHHETDDVVAAKVQAAYRHGLTPILCCGEGLEVRKEGSQVEHVVAQLRAALDGVTREQAASIVIAYEPIWAIGTGEVATPDDAQEVCAAIRTLLAELYSGDLADGVRILYGGSVKAANVAAIMAQEDVDGALVGGASIDPAEFASICRYRDHLTAG",
    "MRTKMIAGNWKMHHRPQEARAFVEELGRVLWARNELYGPLKEGVAEAVLFPTALSLAAVQDALGDLPVSLGAQNAHWEDHGAFTGEIGAPMLADFGCAYILIGHSERRHLFHETEVELARKLRAVLSTSARCLFCVGELLEEREAGKTHQVLERQLLGALEGVTIPDLTDRFAVAYEPVWAIGTGKTASDGDAEEGCGYIRHLVADRYGQETAQHLQVLYGGSVKPGNTAGLMVQGDIDGLLVGGASLEVPSFVGILEAALGILRP",
  ];
  // const aminos = [
  //   "MKPIIAVNYKAYYPYSFGENALRIARDAKRVWEETGVEVILAPPFTEIYRVLKEVEGSGVKVFAQHAVMKDMALNALKAVS",
  //   "MAPRKFFVGGNWKMNGDKKSLGELIHTLNGAKLSADTEVVCGAPSIYLDFARQKLDAKIGVAAQNCYCKELASQHDVDGFLVGGASLKPEFVDIINAKH",
  //   "MRRYLIAGNWKMNTSLETGTALASGLADHVRGRDLPVDVLVCPPFPYLAAVKATAGEAGISVGAQNCYFEASGAFTGEVSVDMLKDIGCDS",
  // ];
  const data = {
    aminos: aminos,
  };

  const getTree = async () => {
    if (aminos.length <= 1) {
      console.log("2個以上のデータを指定してください");
      return;
    }
    const res = await fetch("/tree", {
      method: "POST",
      headers: {
        "Content-Type": "application/json",
      },
      body: JSON.stringify(data),
    });
    const tree = await res.json();
    let nowId = 0;

    const drawTree = (tree, id, parentScore) => {
      const y0 = (0.5 - parentScore) * HEIGHT * 2;
      const y1 = (0.5 - tree[id].score) * HEIGHT * 2;
      // 葉の場合
      if (tree[id].left === -1 && tree[id].right === -1) {
        const x =
          (nowId * (WIDTH - PADDING)) / (aminos.length - 1) + PADDING / 2;
        context.fillRect(x, y0, LINE_WIDTH, y1 - y0);
        nowId++;
        return x;
      }
      // 節
      const leftx = drawTree(tree, tree[id].left, tree[id].score);
      const rightx = drawTree(tree, tree[id].right, tree[id].score);
      const x = (leftx + rightx) / 2;
      context.fillRect(x, y0, LINE_WIDTH, y1 - y0);
      context.fillRect(leftx, y1, rightx - leftx, LINE_WIDTH);
      return x;
    };
    drawTree(tree, tree.length - 1, 1);
  };

  getTree();
};
