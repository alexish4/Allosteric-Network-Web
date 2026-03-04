import React from 'react';
import { Link } from 'react-router-dom';
import './App.css';

const pages = [
  {
    title: "Current Flow Allostery",
    path: "/newallosteric",
    img: "AllostericThumbnail.png"
  },
  {
    title: "PDB Compare",
    path: "/pdbpaircompare",
    img: "PDBCompareThumbnail.png"
  },
  {
    title: "CholBindNet",
    path: "/cholbindnet",
    img: "CholBindThumbnail.png"
  }
];

function Home() {
  return (
    <div className="home-container">
      <h1>Comp Bio Tools</h1>
      <h2>Select a Page</h2>

      <div className="card-grid">
        {pages.map((page, index) => (
          <Link to={page.path} key={index} className="card">
            <img src={page.img} alt={page.title} />
            <div className="card-title">{page.title}</div>
          </Link>
        ))}
      </div>
    </div>
  );
}

export default Home;
