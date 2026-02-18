// Import necessary components
import React from 'react';
import { BrowserRouter as Router, Route, Routes, Link, useLocation } from 'react-router-dom';
import Home from './Home';
import Allosteric from './Allosteric'; 
import PDBCompare from './PDBCompare';
import './App.css';
import NewAllosteric from './NewAllosteric';
import CholNet from './CholNet';  

function App() {
  return (
    <Router>
      <Routes>
        <Route path="/" element={<Home />} />
        <Route path="/allosteric" element={<Allosteric />} />
        <Route path="/pdbpaircompare" element={<PDBCompare />} />
        <Route path="/newallosteric" element={<NewAllosteric />} />
        <Route path="/cholbindnet" element={<CholNet />} />
      </Routes>
      <NavBar />
    </Router>
  );
}

function NavBar() {
  const location = useLocation(); // Get the current URL

  return (
    <nav>
      {location.pathname !== "/" && <Link to="/">Go Back to Home Page |</Link>}
      {location.pathname !== "/newallosteric" && <Link to="/newallosteric"> Allosteric Page</Link>}
      {location.pathname !== "/pdbpaircompare" && <Link to="/pdbpaircompare"> | PDB Compare Page</Link>} 
      {location.pathname !== "/cholbindnet" && <Link to="/cholbindnet"> | CholBindNet Page</Link>}
    </nav>
  );
}

export default App;